source("code/00_functions.R")
library(readr)
library(tidyverse)
library(janitor)
library(ISOweek)
# source("R/00_functions.R")

# files downloaded from Datasus, the source is the Ministry of Health
# data 2015-2020
# https://opendatasus.saude.gov.br/dataset/sim-1979-2019
# data 2021
# https://dados.gov.br/dataset/sistema-de-informacao-sobre-mortalidade

links <- paste0("data_input/brazil/Mortalidade_Geral_", 2000:2020, ".csv")
i <- links[1]
out <- list()
for (i in links){
  
  cat(i)
  # test <- out[[i]]
  # unique(test$SEXO)
  out[[i]] <- 
    read_delim(i, 
               delim = ";",
               col_types = cols(.default = "c")) %>%  
    filter(TIPOBITO == "2") %>% 
    select(date_e = DTOBITO, 
           IDADE,
           SEXO,
           CAUSABAS) %>% 
    mutate(date_e = dmy(date_e),
           year = year(date_e),
           month = month(date_e),
           sex = recode(SEXO,
                        "0" = "UNK",
                        "1" = "m",
                        "2" = "f",
                        "9" = "UNK"),
           age = case_when(IDADE <= 400 ~ "0",
                           IDADE > 400 & IDADE < 500 ~ str_sub(IDADE, 2, 3),
                           IDADE >= 500 & IDADE <= 600 ~ "100",
                           TRUE ~ "UNK"),
           cause_let = str_sub(CAUSABAS, 1, 1),
           cause_num = str_sub(CAUSABAS, 2, 3) %>% as.integer(),
           cause = case_when(
             cause_let == "J" & cause_num %in% 10:18 ~ "pi", 
             cause_let == "J" & !cause_num %in% 10:18 ~ "res",
             cause_let == "I" ~ "cvd",
             TRUE ~ "oth")) %>% 
    group_by(year, month, cause, sex, age) %>% 
    summarize(dts = n(), .groups = "drop")
}

dts <- 
  out %>% 
  bind_rows() %>% 
  ungroup() 

# imputing unknown sex
# ~~~~~~~~~~~~~~~~~~~~
dts2 <- 
  dts %>% 
  summarise(dts = sum(dts), .by = c(year, month, cause, sex, age)) %>% 
  spread(sex, dts) %>% 
  replace_na(list(UNK = 0,
                  m = 0,
                  f = 0)) %>% 
  mutate(t = f + m + UNK,
         t_s = f + m,
         f = (t/t_s)*f,
         m = (t/t_s)*m) %>% 
  select(-UNK, -t_s) %>% 
  gather(f, m, t, key = sex, value = dts)

unique(dts2$age)

# imputing unknown age
# ~~~~~~~~~~~~~~~~~~~~
dts3 <- 
  dts2 %>% 
  drop_na(year) %>% 
  mutate(dts_tot = sum(dts), .by = c(year, month, cause, sex)) %>% 
  filter(age != "UNK") %>% 
  mutate(dts_sum = sum(dts), 
         dts = dts*(dts_tot/dts_sum), 
         .by = c(year, month, cause, sex)) %>% 
  select(-dts_tot, -dts_sum) %>% 
  pivot_wider(names_from = age, values_from = dts,
              values_fill = 0) %>% 
  pivot_longer(!c(year, month, cause, sex), 
               names_to = "age", values_to = "dts") %>% 
  mutate(age = age %>% as.double()) %>% 
  arrange(year, month, cause, sex, age)

# adding all-cause mortality
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
dts4 <- 
  dts3 %>% 
  bind_rows(
    dts3 %>% 
      summarise(dts = sum(dts), .by = c(year, month, sex, age)) %>% 
      mutate(cause = "total")
  ) %>% 
  mutate(date = make_date(d = 15, m = month, y = year))



write_rds(dts4, "data_inter/brazil_monthly_cause_death.rds", compress = "xz")
dt <- read_rds("data_inter/brazil_monthly_cause_death.rds")
# rm(-dt)
unique(dt$cause)
