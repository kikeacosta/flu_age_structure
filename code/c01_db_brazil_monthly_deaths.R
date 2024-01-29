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

# dt %>% 
#   filter(age %in% 50:54,
#          sex == "t",
#          cause == "pi",
#          year %in% 2000:2019) %>% 
#   ggplot()+
#   geom_line(aes(date, dts, col = age, group = age))+
#   theme_bw()
# 
# dt %>% 
#   filter(age %in% 50:54,
#          sex == "t",
#          cause == "total",
#          year %in% 2000:2019) %>% 
#   ggplot()+
#   geom_line(aes(date, dts, col = age, group = age))+
#   theme_bw()
# 
# 
# dt_sml <- 
#   dt %>% 
#   filter(sex == "t") %>% 
#   mutate(age = age - age %% 5) %>% 
#   summarise(dts = sum(dts), .by = c(year, month, date, cause, sex, age))
# 
# dt_sml %>% 
#   filter(age %in% 30:70,
#          sex == "t",
#          # cause %in% c("pi", "total"),
#          year %in% 2000:2019) %>% 
#   ggplot()+
#   geom_line(aes(date, dts, col = age, group = age))+
#   facet_wrap(~cause, scales = "free_y")+
#   theme_bw()
# 
# dt_sml %>% 
#   filter(age %in% 30:70,
#          sex == "t",
#          # cause %in% c("pi", "total"),
#          year %in% 2006:2014) %>% 
#   ggplot()+
#   geom_line(aes(date, dts, col = age, group = age))+
#   facet_wrap(~cause, scales = "free_y")+
#   theme_bw()
# 
# dt_sml %>% 
#   filter(age == 70,
#          sex == "t",
#          # cause %in% c("pi", "total"),
#          year %in% 2000:2019) %>% 
#   ggplot()+
#   geom_line(aes(month, dts, col = year, group = year))+
#   scale_x_continuous(breaks = 1:12)+
#   facet_wrap(~cause, scales = "free_y")+
#   theme_bw()
# 
# 
# chunk <- 
#   dt %>% 
#   filter(sex == "t",
#          cause == "total",
#          age == 60) %>% 
#   mutate(t = 1:n(), 
#          w = ifelse(month %in% 5:9, 0, 1),
#          .by = c(cause, sex, age))
# 
# # 
# 
# mod <- 
#   gam(dts ~ t + s(month) + offset(pop))
# 
# 
# 
# 
# dts3 <- 
#   dts2 %>% 
#   filter(year %in% 2004:2016)
# 
# 
# tot_age <- 
#   dts %>% 
#   group_by(year, week, cause, sex) %>% 
#   summarise(dts = sum(dts), .groups = "drop") %>% 
#   mutate(age = "TOT")
# 
# # re-scaling age and sex
# dts2 <- 
#   dts %>% 
#   filter(Age != "UNK") %>% 
#   bind_rows(tot_age) %>% 
#   group_by(Sex, Year) %>% 
#   do(rescale_age(chunk = .data)) %>% 
#   ungroup() %>%
#   group_by(Age, Year) %>%
#   do(rescale_sex(chunk = .data)) %>%
#   ungroup()
# 
# dts3 <- 
#   dts2 %>% 
#   mutate(Age = Age %>% as.double(),
#          Country = "Brazil",
#          Code = "BRA",
#          Source = "brazil_sim") %>% 
#   select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
#   arrange(Year, Sex, Age)
# 
# # write_rds(dts3, "data_inter/brazil.rds")
# # dts3 <- read_rds("Output/brazil.rds")
# write_csv(dts3, "data_inter/brazil.csv")
# sub