rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)


dt_in <- NULL
dt_in[1] <- 
  read_delim("data_input/brazil/INFLUD09.csv", delim = ";")

for(i in 2:11){
  dt_in[i] <- 
    read_delim(paste0("data_input/brazil/INFLUD", i + 8, ".csv"), delim = ";")
}

dt_in[[1]]

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

in09_11 <- 
  bind_rows(in09, in10, in11)

in09_112 <- 
  in09_11 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, 
         PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         EVOLUCAO, DT_OBITO)

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
         pcr_test = ifelse(PCR_ETIOL %in% 1:4, 1, 0),
         pcr_test2 = ifelse(is.na(pcr_test), 0, 1),
         pcr_test2 = ifelse((pcr_test == 1 | typ != "other" | sub != "other"), 1, 0),
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
         year = year(date_flu)) %>% 
  filter(pcr_test2 == 1) %>% 
  select(year, sex, date_bth, date_flu, typ, sub, outcome, date_dth, date_not,
         pcr_test, pcr_test2)


flu0911 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, typ) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, outcome) %>% 
  summarise(n = n())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

in122 <- 
  in12 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, 
         PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         EVOLUCAO, DT_OBITO)

flu12 <- 
  in122 %>% 
  mutate(sub = case_when(PCR_TIPO_H == 1 ~ "h1",
                         PCR_TIPO_H == 3 ~ "h3",
                         PCR_ETIOL == 3 ~ "b",
                         TRUE ~ "other"),
         sub2 = case_when(RES_FLUASU == 1 ~ "h1",
                          RES_FLUASU == 2 ~ "h1pre",
                          RES_FLUASU == 3 ~ "h3",
                          RES_FLUASU == 4 ~ "an",
                          RES_FLUASU == 5 ~ "h3v",
                          RES_FLUASU == 6 ~ "ao",
                          RES_FLUB == 1 ~ "b",
                          TRUE ~ "other"),
         typ = case_when(RES_FLUA == 1 | 
                           sub %in% c("h1", "h3") |
                           sub %in% c("h1", "h1pre", "h3", "an", "h3v", "ao")
                         ~ "a",
                         RES_FLUB == 1 | 
                           sub %in% c("b") |
                           sub %in% c("b")
                         ~ "b",
                         TRUE ~ "other"),
         pcr_test = ifelse(PCR_ETIOL %in% 1:4 | 
                             RES_FLUA == 1 | 
                             RES_FLUB == 1 |
                             sub != "other" | 
                             sub2 != "other" | 
                             typ %in% c("a", "b"), 
                             1, 0),
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
         year = year(date_flu)) %>% 
  filter(pcr_test == 1) %>% 
  select(year, sex, date_bth, date_flu, typ, sub, sub2, outcome, date_dth, date_not,
         pcr_test)

flu12 %>% 
  group_by(sub) %>% 
  summarise(n = n())

flu12 %>% 
  group_by(sub2) %>% 
  summarise(n = n())

flu12 %>% 
  group_by(typ) %>% 
  summarise(n = n())

flu12 %>% 
  group_by(outcome) %>% 
  summarise(n = n())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

in1318 <- 
  bind_rows(in13, in14, in15, in16, in17, in18) %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, 
         HOSPITAL, 
         RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB,
         EVOLUCAO, DT_OBITO)

flu1318 <- 
  in1318 %>% 
  filter(RES_FLUA == 1 | RES_FLUB == 1) %>% 
  mutate(typ = case_when(RES_FLUA == 1 ~ "a",
                         RES_FLUB == 1 ~ "b",
                         TRUE ~ "other"),
         sub = case_when(RES_FLUASU == 1 ~ "h1",
                         RES_FLUASU == 2 ~ "h1pre",
                         RES_FLUASU == 3 ~ "h3",
                         RES_FLUASU == 4 ~ "an",
                         RES_FLUASU == 5 ~ "h3v",
                         RES_FLUASU == 6 ~ "ao",
                         RES_FLUB == 1 ~ "b"),
       date_bth = dmy(DT_NASC),
       date_flu = dmy(DT_SIN_PRI),
       outcome = case_when(EVOLUCAO == 1 ~ "survived",
                           EVOLUCAO == 2 ~ "death_flu",
                           EVOLUCAO == 3 ~ "death_oth",
                           EVOLUCAO == 4 ~ "death_inv",
                           EVOLUCAO == 5 ~ "missing",
                           TRUE ~ "oth"),
       date_dth = dmy(DT_OBITO),
       sex = CS_SEXO %>% str_to_lower(),
       year = year(date_flu)) %>% 
  select(year, sex, date_bth, date_flu, typ, sub, outcome, date_dth)

flu1318 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu1318 %>% 
  group_by(year, typ) %>% 
  summarise(n = n())

flu1318 %>% 
  group_by(year, outcome) %>% 
  summarise(n = n())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


in192 <- 
  in19 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, 
         POS_PCRFLU,POS_IF_FLU, TP_FLU_IF,
         TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,
         EVOLUCAO, DT_EVOLUCA)

flu19 <- 
  in192 %>% 
  filter(POS_PCRFLU == 1) %>% 
  mutate(typ = case_when(TP_FLU_PCR == 1 ~ "a",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         sub = case_when(PCR_FLUASU == 1 ~ "h1",
                         PCR_FLUASU == 2 ~ "h3",
                         PCR_FLUASU %in% 3:6 ~ "an",
                         FLUASU_OUT == 1 ~ "ao",
                         TP_FLU_PCR == 2 ~ "b"),
         date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         date_dth = dmy(DT_EVOLUCA),
         sex = CS_SEXO %>% str_to_lower(),
         year = year(date_flu)) %>% 
  select(year, sex, date_bth, date_flu, typ, sub, outcome, date_dth)

flu19 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(outcome) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(EVOLUCAO) %>% 
  summarise(n = n())

