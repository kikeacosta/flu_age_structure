rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

# Brazil
# ~~~~~~

# to look at the cohort composition of each flu status
library(wpp2022)
data(popAge1dt)
br_pop <- 
  popAge1dt %>% 
  filter(name == "Brazil",
         year %in% c(2010, 2016, 2019)) %>% 
  mutate(cohort = year - age) %>% 
  select(year, cohort, n = pop) %>% 
  group_by(year) %>% 
  mutate(cx = n/sum(n),
         status = "all brazil wpp")


# data on flu survellance
# Banco de Dados de Síndrome Respiratória Aguda Grave
# Notificações de Síndrome Gripal
# https://opendatasus.saude.gov.br/organization/ministerio-da-saude
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
in09 <- read_delim("data_input/brazil/INFLUD09.csv", delim = ";")
in10 <- read_delim("data_input/brazil/INFLUD10.csv", delim = ";")
in11 <- read_delim("data_input/brazil/INFLUD11.csv", delim = ";")
in12 <- read_delim("data_input/brazil/INFLUD12.csv", delim = ";")
in13 <- read_delim("data_input/brazil/INFLUD13.csv", delim = ";")
in14 <- read_delim("data_input/brazil/INFLUD14.csv", delim = ";")
in15 <- read_delim("data_input/brazil/INFLUD15.csv", delim = ";")
in16 <- read_delim("data_input/brazil/INFLUD16.csv", delim = ";")
in17 <- read_delim("data_input/brazil/INFLUD17.csv", delim = ";")
in18 <- read_delim("data_input/brazil/INFLUD18.csv", delim = ";")
in19 <- read_delim("data_input/brazil/INFLUD19.csv", delim = ";")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unique(in09$NU_IDADE_N) %>% sort()
unique(in11$NU_IDADE_N) %>% sort()
unique(in18$NU_IDADE_N) %>% sort()
unique(in19$NU_IDADE_N) %>% sort()

unique(in09$PCR_ETIOL) %>% sort

test <- 
  in09 %>% 
  select(1, PCR_RES, PCR_ETIOL, PCR_TIPO_H, PCR_TIPO_N, RES_FLUA, RES_FLUASU,
         DS_OUTSUB)

test2 <- 
  test %>% 
  filter(PCR_ETIOL %in% 1:4)

test3 <- 
  test %>% 
  filter(PCR_RES == 1,
         PCR_ETIOL == 5)

test4 <- 
  test %>% 
  filter(is.na(PCR_RES))

# ~~~~~~~~~~~~~~~~~~~
# data 2009-2011 ====
# ~~~~~~~~~~~~~~~~~~~

in09_11 <- 
  bind_rows(in09, in10, in11)

colnames(in09) %>% sort()
colnames(in10) %>% sort()
colnames(in11) %>% sort()
colnames(in12) %>% sort()
colnames(in13) %>% sort()
colnames(in18) %>% sort()

in09_112 <- 
  in09_11 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, NU_IDADE_N,
         HOSPITAL, 
         PCR_RES,
         PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         VACINA, 
         EVOLUCAO, DT_OBITO)

unique(in09_112$NU_IDADE_N) %>% sort()

in09_113 <- 
  in09_112 %>% 
  mutate(pcr = case_when(PCR_RES == 1 ~ "pos",
                          PCR_RES == 2 ~ "neg",
                          PCR_RES == 3 ~ "inc",
                          PCR_RES == 4 ~ "not",
                          TRUE ~ NA_character_),
         typ = case_when(PCR_ETIOL == 1 ~ "a",
                         PCR_ETIOL == 2 ~ "a",
                         PCR_ETIOL == 3 ~ "b",
                         PCR_ETIOL == 4 ~ "av",
                         PCR_ETIOL == 5 ~ "other",
                         TRUE ~ NA_character_),
         sub = case_when(PCR_TIPO_H == 1 ~ "h1",
                         PCR_TIPO_H == 3 ~ "h3",
                         PCR_ETIOL == 3 ~ "b",
                         TRUE ~ "other"),
         typ = case_when(!is.na(PCR_TIPO_H) | !is.na(PCR_TIPO_N) ~ "a",
                          sub %in% c("b") ~ "b",
                          TRUE ~ typ),
         flu = ifelse(PCR_RES == 1 | 
                           PCR_ETIOL %in% 1:4 | 
                           !is.na(PCR_TIPO_H) | 
                           !is.na(PCR_TIPO_N), 1, 0), 
         date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         date_dth = dmy(DT_OBITO),
         sex = CS_SEXO %>% str_to_lower(),
         age = case_when(NU_IDADE_N < 4000 ~ 0,
                         NU_IDADE_N > 4000 ~ NU_IDADE_N - 4000),
         year = year(date_flu),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk"),
         state = SG_UF_NOT) %>% 
  select(year, state, sex, 
         age, date_bth, flu, 
         date_flu, typ, sub, 
         outcome, date_dth, 
         date_not, vaccine,
         PCR_RES, PCR_ETIOL, 
         PCR_TIPO_H, PCR_TIPO_N)

all <- 
  in09_113 %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "all SRAG")

miss <- 
  in09_113 %>% 
  filter(is.na(flu)) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "nas")

flus <- 
  in09_113 %>% 
  filter(flu == 1) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "flu")

bind_rows(all, miss, flus, br_pop %>% filter(year == 2010) %>% select(cohort, cx, status)) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()

ggsave("figures/explor_data/cohort_composition_09_11.png",
       w = 8, h = 3)

bind_rows(all, miss, flus, br_pop %>% filter(year == 2010) %>% select(cohort, cx, status)) %>% 
  filter(cohort %in% 1930:2000) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  theme_bw()
ggsave("figures/explor_data/cohort_composition_09_11_zoom.png",
       w = 8, h = 3)


flu0911 <- 
  in09_113 %>% 
  filter(flu == 1) %>% 
  select(year, state, sex, age, date_bth, date_flu, typ, sub, 
         outcome, date_dth, 
         date_not, vaccine)


flu0911 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, typ) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, outcome) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, vaccine) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, age) %>% 
  summarise(n = n())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~
# data 2012-2018 ====
# ~~~~~~~~~~~~~~~~~~~

in12_18 <- 
  bind_rows(in12, in13, in14, in15, in16, in17, in18) %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, NU_IDADE_N,
         HOSPITAL, 
         PCR_RES, PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         VACINA, DT_UT_DOSE,
         EVOLUCAO, DT_OBITO)


in12_182 <- 
  in12_18 %>% 
  mutate(test1 = case_when(PCR_RES == 1 ~ "pos",
                         PCR_RES == 2 ~ "neg",
                         PCR_RES == 3 ~ "inc",
                         PCR_RES == 4 ~ "not",
                         TRUE ~ NA_character_),
         test2 = case_when(RES_FLUA == 1 | RES_FLUB == 1 ~ "pos",
                          RES_FLUA == 2 & RES_FLUB == 2 ~ "neg",
                          RES_FLUA == 3 & RES_FLUB == 3 ~ "inc",
                          RES_FLUA == 4 & RES_FLUB == 4 ~ "not",
                          TRUE ~ NA_character_),
         test = case_when(test1 == "pos" | test2 == "pos" ~ "pos",
                          test1 == "neg" & test2 == "neg" ~ "neg",
                          test1 == "inc" & test2 == "inc" ~ "inc",
                          test1 == "not" & test2 == "not" ~ "not",
                          TRUE ~ NA_character_),
         sub1 = case_when(PCR_TIPO_H == 1 ~ "h1",
                          PCR_TIPO_H == 3 ~ "h3",
                          is.na(PCR_TIPO_H) ~ NA_character_,
                          TRUE ~ "other"),
         sub2 = case_when(RES_FLUASU == 1 ~ "h1",
                          RES_FLUASU == 2 ~ "h1pre",
                          RES_FLUASU == 3 ~ "h3",
                          RES_FLUASU == 4 ~ "an",
                          RES_FLUASU == 5 ~ "h3v",
                          RES_FLUASU == 6 ~ "ao",
                          RES_FLUB == 1 ~ "b",
                          TRUE ~ NA_character_),
         sub = case_when((!is.na(sub1) & sub1 != "other") & (is.na(sub2) | sub2 == "other") ~ sub1,
                         (is.na(sub1) & sub1 == "other") & (!is.na(sub2) | sub2 != "other") ~ sub2,
                         (sub1 == "h1" | sub2 == "h1") & (sub2 != "h1pre") ~ "h1",
                         sub2 == "h1pre" ~ "h1pre",
                         (sub1 == "h3" | sub2 == "h3") & (sub2 != "h3v") ~ "h3",
                         sub2 == "h3v" ~ "h3v",
                         sub1 == "b" | sub2 == "b" ~ "b",
                         is.na(sub1) & is.na(sub2) ~ NA_character_,
                         TRUE ~ sub2),
         typ1 = case_when(PCR_ETIOL == 1 ~ "a",
                          PCR_ETIOL == 2 ~ "a",
                          PCR_ETIOL == 3 ~ "b",
                          PCR_ETIOL == 4 ~ "av",
                          PCR_ETIOL == 5 ~ "other",
                          TRUE ~ NA_character_),
         typ2 = case_when(RES_FLUA == 1 ~ "a",
                          RES_FLUB == 1 ~ "b",
                          TRUE ~ "other"),
         typ = case_when(typ1 == "a" | typ2 == "a" ~ "a",
                         typ1 == "b" | typ2 == "b" ~ "b",
                         sub1 %in% c("h1", "h3") |
                           sub2 %in% c("h1", "h1pre", "h3", "an", "h3v", "ao") ~ "a",
                         sub1 %in% c("b") |
                           sub2 %in% c("b") ~ "b",
                         TRUE ~ typ1)) %>% 
  mutate(sub = case_when(typ == "b" & sub == "other" ~ "b", 
                         typ == "a" & sub == "other" ~ "an", 
                         typ %in% c("new", "av") & sub == "other" ~ typ,
                         TRUE ~ sub),
         flu = ifelse(PCR_RES == 1 | 
                        PCR_ETIOL %in% 1:4 | 
                        !is.na(PCR_TIPO_H) | 
                        !is.na(PCR_TIPO_N) |
                        RES_FLUA == 1 | 
                        RES_FLUB == 1 | 
                        typ %in% c("a", "b"), 
                      1, 0), 
         date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         date_vac = dmy(DT_UT_DOSE),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         date_dth = dmy(DT_OBITO),
         sex = CS_SEXO %>% str_to_lower(),
         age = case_when(NU_IDADE_N < 4000 ~ 0,
                         NU_IDADE_N > 4000 ~ NU_IDADE_N - 4000),
         year = year(date_flu),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk")) 

all12_18 <- 
  in12_182 %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "all SRAG")

miss12_18 <- 
  in12_182 %>% 
  filter(is.na(flu)) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "nas")

misssub12_18 <- 
  in12_182 %>% 
  filter(is.na(sub)) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "missup")

flus12_18 <- 
  in12_182 %>% 
  filter(flu == 1) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "flu")

bind_rows(all12_18, 
          miss12_18, 
          flus12_18, 
          br_pop %>% 
            filter(year == 2016) %>% 
            select(cohort, cx, status)) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()

ggsave("figures/explor_data/cohort_composition_12_18.png",
       w = 8, h = 3)

bind_rows(all12_18, 
          miss12_18, 
          # misssub12_18,
          flus12_18, 
          br_pop %>% 
            filter(year == 2016) %>% 
            select(cohort, cx, status)) %>% 
  filter(cohort %in% 1930:2000) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  theme_bw()
ggsave("figures/explor_data/cohort_composition_12_18_zoom.png",
       w = 8, h = 3)


flu1218 <- 
  in12_182 %>% 
  filter(flu == 1) %>% 
  select(year, sex, age, date_bth, date_flu, typ1, typ2, typ, sub1, sub2, sub, 
         outcome, date_dth, date_not, vaccine, date_vac) %>% 
  select(-sub1, -sub2, -typ1, -typ2)

flu1218 %>% 
  group_by(year, sub) %>% 
  summarise(n = n()) %>% 
  spread(sub, n)

flu1218 %>% 
  group_by(year, typ) %>% 
  summarise(n = n()) %>% 
  spread(typ, n)

flu1218 %>% 
  group_by(year, typ, sub) %>% 
  summarise(n = n()) 

flu1218 %>% 
  group_by(year, age) %>% 
  summarise(n = n()) 



# ~~~~~~~~~~~~~~
# data 2019 ====
# ~~~~~~~~~~~~~~

in192 <- 
  in19 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, NU_IDADE_N,
         HOSPITAL, 
         PCR_RESUL, POS_PCRFLU,POS_IF_FLU, TP_FLU_IF,
         TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,
         VACINA, DT_UT_DOSE,
         EVOLUCAO, DT_EVOLUCA)

in193 <- 
  in192 %>% 
  mutate(test = case_when(POS_PCRFLU == 1 ~ "pos",
                          POS_PCRFLU == 2 ~ "neg",
                          POS_PCRFLU == 3 ~ "inc",
                          PCR_RESUL == 4 ~ "not"),
         typ = case_when(TP_FLU_PCR == 1 ~ "a",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ NA_character_),
         sub = case_when(PCR_FLUASU == 1 ~ "h1",
                         PCR_FLUASU == 2 ~ "h3",
                         PCR_FLUASU %in% 3:5 ~ "an",
                         PCR_FLUASU == 6 | FLUASU_OUT == 1 ~ "ao",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ NA_character_),
         date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         date_vac = dmy(DT_UT_DOSE),
         flu = ifelse(POS_PCRFLU == 1 | 
                             typ != "other" | 
                             sub != "other", 
                           1, 0),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         date_dth = dmy(DT_EVOLUCA),
         sex = CS_SEXO %>% str_to_lower(),
         age = NU_IDADE_N,
         year = year(date_flu),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk"))


all19 <- 
  in193 %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "all SRAG")

miss19 <- 
  in193 %>% 
  filter(is.na(flu)) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "nas")

misssub19 <- 
  in193 %>% 
  filter(is.na(sub)) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "missup")

flus19 <- 
  in193 %>% 
  filter(flu == 1) %>% 
  mutate(cohort = year(date_bth)) %>% 
  group_by(cohort) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(cx = n / sum(n), 
         status = "flu")

bind_rows(all19, 
          miss19, 
          flus19, 
          br_pop %>% 
            filter(year == 2019) %>% 
            select(cohort, cx, status)) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()

ggsave("figures/explor_data/cohort_composition_19.png",
       w = 8, h = 3)

bind_rows(all19, 
          miss19, 
          flus19, 
          br_pop %>% 
            filter(year == 2019) %>% 
            select(cohort, cx, status)) %>% 
  filter(cohort %in% 1930:2000) %>% 
  ggplot()+
  geom_line(aes(cohort, cx, col = status))+
  theme_bw()
ggsave("figures/explor_data/cohort_composition_19_zoom.png",
       w = 8, h = 3)




flu19 <- 
  in193 %>% 
  filter(flu == 1) %>% 
  select(year, sex, age, date_bth, date_flu, date_not, typ, sub, outcome, date_dth, vaccine, date_vac)

flu19 %>% 
  group_by(year, typ) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(outcome) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(vaccine) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(age) %>% 
  summarise(n = n())


# ~~~~~~~~~~~~~~~~~~~~~~
# All data together ====
# ~~~~~~~~~~~~~~~~~~~~~~
flu <- 
  bind_rows(flu0911, flu1218, flu19)

test <- 
  flu %>% 
  group_by(year, age) %>% 
  summarise(n = n())


write_rds(flu, "data_inter/flu_data_brazil_2009_2019_v2.rds")

# comparison with previous version
flu_prev <- read_rds("data_inter/flu_data_brazil_2009_2019.rds")
# flu <- read_rds("data_inter/flu_data_brazil_2009_2019.rds")

sub_prev <- 
  flu_prev %>% 
  group_by(sub) %>% 
  summarise(n = n()) %>%  
  ungroup()

sub_curr <- 
  flu %>% 
  group_by(sub) %>% 
  summarise(n = n()) %>%  
  ungroup()




















vac_date <- 
  flu %>% 
  filter(year %in% 2013:2019) %>% 
  mutate(test = ifelse(is.na(date_vac), "n", "y")) %>% 
  group_by(year, test) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  spread(test, n)

# table(vac_date$is_vac_date)



# 
# 
# # 09-11
# select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
#        DT_NASC, CS_SEXO, NU_IDADE_N,
#        HOSPITAL, 
#        PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
#        PCR_TIPO_H, PCR_TIPO_N,
#        VACINA, 
#        EVOLUCAO, DT_OBITO)
# 
# # 12
# select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
#        DT_NASC, CS_SEXO, NU_IDADE_N,
#        HOSPITAL, 
#        PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
#        PCR_TIPO_H, PCR_TIPO_N,
#        VACINA, DT_UT_DOSE,
#        EVOLUCAO, DT_OBITO)
# 
# # 13-18
# select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
#        DT_NASC, CS_SEXO, NU_IDADE_N,
#        HOSPITAL, 
#        PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
#        PCR_TIPO_H, PCR_TIPO_N,
#        VACINA, DT_UT_DOSE,
#        EVOLUCAO, DT_OBITO)
# 
# # 19
# select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
#        DT_NASC, CS_SEXO, NU_IDADE_N,
#        HOSPITAL, 
#        POS_PCRFLU,POS_IF_FLU, TP_FLU_IF,
#        TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,
#        VACINA, DT_UT_DOSE,
#        EVOLUCAO, DT_EVOLUCA)
