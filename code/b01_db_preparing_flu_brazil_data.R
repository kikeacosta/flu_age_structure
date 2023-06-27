rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

# Brazil

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

unique(in11$NU_IDADE_N) %>% sort()
unique(in18$NU_IDADE_N) %>% sort()
unique(in19$NU_IDADE_N) %>% sort()


test <- 
  in12 %>% 
  select(1, PCR_ETIOL, RES_FLUA, RES_FLUASU)


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
         PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         VACINA, 
         EVOLUCAO, DT_OBITO)

unique(in09_112$NU_IDADE_N) %>% sort()

flu0911 <- 
  in09_112 %>% 
  mutate(typ = case_when(PCR_ETIOL == 1 ~ "new",
                         PCR_ETIOL == 2 ~ "a",
                         PCR_ETIOL == 3 ~ "b",
                         PCR_ETIOL == 4 ~ "av",
                         TRUE ~ "other"),
         sub = case_when(PCR_TIPO_H == 1 ~ "h1",
                         PCR_TIPO_H == 3 ~ "h3",
                         PCR_ETIOL == 3 ~ "b",
                         TRUE ~ "other"),
         typ = case_when(sub %in% c("h1", "h3") ~ "a",
                         sub %in% c("b") ~ "b",
                         TRUE ~ typ),
         pcr_test = ifelse(!is.na(PCR_ETIOL) & PCR_ETIOL %in% 1:4, 1, 0),
         flu = ifelse((pcr_test == 1 | typ != "other" | sub != "other"), 1, 0),
         date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         # date_vac = dmy(DT_UT_DOSE),
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

# ~~~~~~~~~~~~~~~~~~~
# data 2012-2018 ====
# ~~~~~~~~~~~~~~~~~~~

in1218 <- 
  bind_rows(in12, in13, in14, in15, in16, in17, in18) %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, NU_IDADE_N,
         HOSPITAL, 
         PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         VACINA, DT_UT_DOSE,
         EVOLUCAO, DT_OBITO)

flu1218 <- 
  in1218 %>% 
  mutate(sub1 = case_when(PCR_TIPO_H == 1 ~ "h1",
                          PCR_TIPO_H == 3 ~ "h3",
                          TRUE ~ "other"),
         sub2 = case_when(RES_FLUASU == 1 ~ "h1",
                          RES_FLUASU == 2 ~ "h1pre",
                          RES_FLUASU == 3 ~ "h3",
                          RES_FLUASU == 4 ~ "an",
                          RES_FLUASU == 5 ~ "h3v",
                          RES_FLUASU == 6 ~ "ao",
                          RES_FLUB == 1 ~ "b",
                          TRUE ~ "other"),
         sub = case_when(sub1 != "other" & sub2 == "other" ~ sub1,
                         sub1 == "other" & sub2 != "other" ~ sub2,
                         (sub1 == "h1" | sub2 == "h1") & (sub2 != "h1pre") ~ "h1",
                         sub2 == "h1pre" ~ "h1pre",
                         (sub1 == "h3" | sub2 == "h3") & (sub2 != "h3v") ~ "h3",
                         sub2 == "h3v" ~ "h3v",
                         sub1 == "b" | sub2 == "b" ~ "b",
                         TRUE ~ sub2),
         typ1 = case_when(PCR_ETIOL == 1 ~ "new",
                          PCR_ETIOL == 2 ~ "a",
                          PCR_ETIOL == 3 ~ "b",
                          PCR_ETIOL == 4 ~ "av",
                          TRUE ~ "other"),
         typ2 = case_when(RES_FLUA == 1 ~ "a",
                          RES_FLUB == 1 ~ "b",
                          TRUE ~ "other"),
         typ = case_when(typ1 == "a" | typ2 == "a" ~ "a",
                         typ1 == "b" | typ2 == "b" ~ "b",
                         sub1 %in% c("h1", "h3") |
                           sub2 %in% c("h1", "h1pre", "h3", "an", "h3v", "ao") ~ "a",
                         sub1 %in% c("b") |
                           sub2 %in% c("b") ~ "b",
                         TRUE ~ typ1),
         sub = case_when(typ == "b" & sub == "other" ~ "b", 
                         typ == "a" & sub == "other" ~ "an", 
                         typ %in% c("new", "av") & sub == "other" ~ typ,
                         TRUE ~ sub),
         flu = ifelse(PCR_ETIOL %in% 1:4 | 
                        RES_FLUA == 1 | 
                        RES_FLUB == 1 |
                        sub1 != "other" | 
                        sub2 != "other" | 
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
                             TRUE ~ "unk")) %>% 
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
         POS_PCRFLU,POS_IF_FLU, TP_FLU_IF,
         TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,
         VACINA, DT_UT_DOSE,
         EVOLUCAO, DT_EVOLUCA)

flu19 <- 
  in192 %>% 
  mutate(typ = case_when(TP_FLU_PCR == 1 ~ "a",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         sub = case_when(PCR_FLUASU == 1 ~ "h1",
                         PCR_FLUASU == 2 ~ "h3",
                         PCR_FLUASU %in% 3:6 ~ "an",
                         FLUASU_OUT == 1 ~ "ao",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
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
                             TRUE ~ "unk")) %>% 
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


write_rds(flu, "data_inter/flu_data_brazil_2009_2019.rds")



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
