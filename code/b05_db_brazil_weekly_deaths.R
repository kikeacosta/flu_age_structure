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
           cause = ifelse(cause_let == "J" & cause_num %in% 10:18, "pi", "oth"),
           week = date2ISOweek(date_e),
           week = str_sub(week, 1, 8)) %>% 
    group_by(year, week, cause, sex, age) %>% 
    summarize(dts = n(), .groups = "drop") %>% 
    spread(sex, dts) %>% 
    replace_na(list(UNK = 0,
                    m = 0,
                    f = 0)) %>% 
    mutate(t = f + m + UNK) %>% 
    select(-UNK) %>% 
    gather(f, m, t, key = sex, value = dts) %>% 
    group_by(year, week, cause, sex, age) %>% 
    summarize(dts = sum(dts), .groups = "drop")
}

dts <- 
  out %>% 
  bind_rows() %>% 
  ungroup() 

write_rds(dts, "data_inter/brazil_weekly_cause_death.rds")
dts <- read_rds("data_inter/brazil_weekly_cause_death.rds")


dts2 <- 
  dts %>% 
  drop_na(year) %>% 
  mutate(isoweek = paste0(week, "-7"),
         date = ISOweek2date(isoweek)) %>% 
  group_by(date, year, cause, sex, age) %>% 
  summarize(dts = sum(dts), .groups = "drop")
  
unique(dts2$cause)

dts2 %>% 
  filter(age == 80,
         sex == "t",
         cause == "pi",
         year %in% 2004:2016) %>% 
  ggplot()+
  geom_line(aes(date, dts))+
  theme_bw()

dts3 <- 
  dts2 %>% 
  filter(year %in% 2004:2016)

chunk <- 
  dts3 %>% 
  filter(sex == "t", )
  
mod <- 
  gam()





tot_age <- 
  dts %>% 
  group_by(year, week, cause, sex) %>% 
  summarise(dts = sum(dts), .groups = "drop") %>% 
  mutate(age = "TOT")

# re-scaling age and sex
dts2 <- 
  dts %>% 
  filter(Age != "UNK") %>% 
  bind_rows(tot_age) %>% 
  group_by(Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>%
  ungroup()

dts3 <- 
  dts2 %>% 
  mutate(Age = Age %>% as.double(),
         Country = "Brazil",
         Code = "BRA",
         Source = "brazil_sim") %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Year, Sex, Age)

# write_rds(dts3, "data_inter/brazil.rds")
# dts3 <- read_rds("Output/brazil.rds")
write_csv(dts3, "data_inter/brazil.csv")
sub