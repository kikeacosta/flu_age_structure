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
# https://opendatasus-saude-gov-br.translate.goog/dataset/srag-2021-a-2023/resource/dd91a114-47a6-4f21-bcd5-86737d4fc734?_x_tr_sl=auto&_x_tr_tl=fr&_x_tr_hl=fr
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
in20 <- read_delim("data_input/brazil/INFLUD20.csv", delim = ";")
in21 <- read_delim("data_input/brazil/INFLUD21.csv", delim = ";")
in22 <- read_delim("data_input/brazil/INFLUD22.csv", delim = ";")
in23 <- read_delim("data_input/brazil/INFLUD23.csv", delim = ";")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

in09_112 <- 
  in09_11 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO, HOSPITAL, 
         PCR_RES,PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU, DS_OUTSUB, 
         PCR_TIPO_H, PCR_TIPO_N,
         HEM_TIPO_H, HEM_TIPO_N,
         VACINA,EVOLUCAO, DT_OBITO, CLASSI_FIN, CLASSI_OUT)

flu0911 <- 
  in09_112 %>% 
  mutate(date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         # date_vac = dmy(DT_UT_DOSE),
         date_not = dmy(DT_NOTIFIC),
         date_dth = dmy(DT_OBITO),
         sex = CS_SEXO %>% str_to_lower(),
         year = year(date_flu),
         pcr_res = case_when (PCR_RES == 1 ~ "pos",
                              PCR_RES == 2 ~ "neg",
                              PCR_RES == 3 ~ "inconclusive",
                              PCR_RES == 4 ~ "notdone",
                              TRUE~"unk"),
         typ = case_when(PCR_ETIOL == 1 ~ "new",
                         PCR_ETIOL == 2 ~ "a",
                         PCR_ETIOL == 3 ~ "b",
                         PCR_ETIOL == 4 ~ "av",
                         TRUE ~ "other"),
         sub = case_when(PCR_ETIOL == 1 ~ "h1",
                         PCR_TIPO_H == 1 ~ "h1",
                         PCR_TIPO_H == 3 ~ "h3",
                         PCR_ETIOL == 3 ~ "b",
                         TRUE ~ "other"),
         classi_fin = case_when(CLASSI_FIN == 1 ~ "influenza",
                                CLASSI_FIN == 2 ~ "other etiol",
                                CLASSI_FIN == 3 ~ "discarded",
                                TRUE ~ "unk"),
         typ = case_when(sub %in% c("h1", "h3") ~ "a",
                         sub %in% c("b") ~ "b",
                         TRUE ~ typ),
         typ = case_when((typ == "other" & classi_fin == "influenza") ~ "ab",
                          TRUE ~ typ),
         sub = case_when((sub != "h1" & HEM_TIPO_H==1) ~ "h1",
                          TRUE ~ sub),
         sub = case_when((sub == "other" & PCR_RES==2) ~ "notflu",
                         TRUE ~ sub),
         sub = case_when((sub == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ sub),
         sub = case_when((pcr_res != "pos" & classi_fin == "influenza" 
                           & month(date_flu)>6 & month(date_flu)<9
                           & year(date_flu)==2009) ~ "new",
                          TRUE ~ sub),
         pcr_test = ifelse(!is.na(PCR_ETIOL) & PCR_ETIOL %in% 1:4, 1, 0),
         flu = ifelse((pcr_test == 1 | typ != "other" | sub != "other"), 1, 0),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk"),
         state = SG_UF_NOT) %>% 
#  filter(flu == 1) %>% 
  select(year, sex, date_bth, date_flu, date_dth, date_not,
         PCR_RES,PCR_ETIOL, RES_FLUA, RES_FLUB, PCR_TIPO_H, PCR_TIPO_N,
         HEM_TIPO_H, HEM_TIPO_N,typ, sub, classi_fin, pcr_res,
         outcome, vaccine)


flu0911_filter<-
  flu0911 %>%
  filter(PCR_RES == 2, PCR_ETIOL == 1 ) %>% 
  arrange(date_bth)

flu0911_filter<-
  flu0911 %>%
  filter(classi_fin == "influenza", pcr_res!="neg")

flu0911_filter<-
  flu0911 %>%
  filter( PCR_RES=="NA")


table0911 <-
  flu0911 %>%
  group_by(PCR_RES, PCR_ETIOL) %>%
  summarise(n = n())

flu0911 %>% 
  group_by(sub) %>% 
  summarise(n = n())


flu0911 %>% 
  group_by(year, outcome) %>% 
  summarise(n = n())

flu0911 %>% 
  group_by(year, vaccine) %>% 
  summarise(n = n())

# ~~~~~~~~~~~~~~~~~~~
# data 2012-2018 ====
# ~~~~~~~~~~~~~~~~~~~

in1218 <- 
  bind_rows(in12, in13, in14, in15, in16, in17, in18) %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, NU_ANO, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, 
         CLASSI_FIN,PCR_RES, PCR_ETIOL, RES_FLUA, RES_FLUB, RES_FLUASU,
         PCR_TIPO_H, PCR_TIPO_N,HEM_TIPO_H, HEM_TIPO_N, 
         DS_OUTSUB, 
         VACINA,
         EVOLUCAO, DT_OBITO)

# table_in1218<-
#   in1218 %>% 
#   group_by(PCR_RES, PCR_ETIOL) %>% 
#   summarise(n = n())


flu1218 <- 
  in1218 %>% 
  mutate(date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         date_dth = dmy(DT_OBITO),
         sex = CS_SEXO %>% str_to_lower(),
         year = year(date_flu),
          
         sub1 = case_when(PCR_TIPO_H == 1 ~ "h1",
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
         typ1 = case_when(PCR_ETIOL == 1 ~ "new", #not in the 2012 manual but in the 2009-11
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
         classi_fin = case_when(CLASSI_FIN == 1 ~ "influenza",
                                CLASSI_FIN == 2 ~ "other_resp",
                                CLASSI_FIN == 3 ~ "other_etiol",
                                CLASSI_FIN == 4 ~ "not_specified",
                                TRUE ~ "unk"),
         flu = ifelse(PCR_ETIOL %in% 1:4 | 
                        RES_FLUA == 1 | 
                        RES_FLUB == 1 |
                        sub1 != "other" | 
                        sub2 != "other" | 
                        typ %in% c("a", "b"), 
                      1, 0),
         # revise type and sub to "ab" when classi_fin==influenza
         typ = case_when((typ == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ typ),
         sub = case_when((sub == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ sub),
         # extract from "other" those with "negative test for influenza"
         typ = case_when((year==2012 & typ == "other" & PCR_RES==2 & PCR_ETIOL!=1) ~ "notflu",
                         TRUE ~ typ),
         sub = case_when((year==2012 & sub == "other" & PCR_RES==2 & PCR_ETIOL!=1) ~ "notflu",
                         TRUE ~ sub),
         typ = case_when((year>=2012 & flu==0) ~ "notflu",
                         TRUE ~ typ),
         sub = case_when((year>=2012 & flu==0) ~ "notflu",
                         TRUE ~ sub),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk")) %>% 
  select(year, sex, date_bth, date_flu,  date_dth, date_not,
         PCR_RES,PCR_ETIOL, RES_FLUA, RES_FLUB, PCR_TIPO_H, PCR_TIPO_N,
         HEM_TIPO_H, HEM_TIPO_N, typ, sub, classi_fin, CLASSI_FIN,
         outcome, vaccine) #%>% 
 # select(-sub1, -sub2, -typ1, -typ2)

flu1218 %>% 
  group_by(sub) %>% 
  summarise(n = n())

flu1218 %>% 
  group_by(year, flu) %>% 
  summarise(n = n())


flu1218 %>% 
  group_by(year, sub) %>% 
  summarise(n = n()) %>% 
  spread(sub, n)

table1218 <-
  flu1218 %>% 
  group_by(PCR_RES, classi_fin) %>% 
  summarise(n = n()) %>% 
  spread(classi_fin, n)

table1218 <-
  flu1218 %>% 
  group_by(PCR_RES, CLASSI_FIN) %>% 
  summarise(n = n()) %>% 
  spread(CLASSI_FIN, n)

flu1218 %>% 
  group_by(year, typ, sub) %>% 
  summarise(n = n()) 


# ~~~~~~~~~~~~~~
# data 2019 ====
# ~~~~~~~~~~~~~~

in19_2 <- 
  in19 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, 
         PCR_RESUL, POS_PCRFLU,POS_IF_FLU, TP_FLU_IF, 
         TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,
         VACINA,
         EVOLUCAO, DT_EVOLUCA, CLASSI_FIN)

flu19 <- 
  in19_2 %>% 
  mutate(date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         date_dth = dmy(DT_EVOLUCA),
         sex = CS_SEXO %>% str_to_lower(),
         year = year(date_flu),
         
         typ = case_when(TP_FLU_PCR == 1 ~ "a",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         sub = case_when(PCR_FLUASU == 1 ~ "h1",
                         PCR_FLUASU == 2 ~ "h3",
                         PCR_FLUASU %in% 3:6 ~ "an",
                         FLUASU_OUT == 1 ~ "ao",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         classi_fin = case_when(CLASSI_FIN == 1 ~ "influenza",
                                CLASSI_FIN == 2 ~ "other_resp",
                                CLASSI_FIN == 3 ~ "other_etiol",
                                CLASSI_FIN == 4 ~ "unk",
                                TRUE ~ "unk"),
         # revise type and sub to "ab" when classi_fin==influenza
         typ = case_when((typ == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ typ),
         sub = case_when((sub == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ sub),
         flu = ifelse(POS_PCRFLU == 1 | 
                             typ != "other" | 
                             sub != "other", 
                           1, 0),
         typ = case_when((typ=="other" & flu==0) ~ "notflu",
                         TRUE ~typ),
         sub = case_when((sub=="other" & flu==0) ~ "notflu",
                         TRUE ~sub),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk")) %>% 
  #filter(flu == 1) %>% 
  select(year, sex, date_bth, date_flu, date_not, date_dth, 
         typ, sub, classi_fin,
         outcome, vaccine)

flu19 %>% 
  group_by(year, sub) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(outcome) %>% 
  summarise(n = n())

flu19 %>% 
  group_by(vaccine) %>% 
  summarise(n = n())


# ~~~~~~~~~~~~~~
# data 2020-23 ====
# ~~~~~~~~~~~~~~

in20_a <-
  in20 %>% 
  mutate(FATOR_RISC = as.double(FATOR_RISC))

in23_a <-
  in23 %>%
  mutate(OBES_IMC = as.character(OBES_IMC))

in2023 <- 
  bind_rows(in20_a,in21,in22,in23_a)

in2023_2 <- 
  in2023 %>% 
  select(SG_UF_NOT, DT_NOTIFIC, SEM_NOT, DT_SIN_PRI, 
         DT_NASC, CS_SEXO,
         HOSPITAL, TP_FLU_AN,
         POS_PCRFLU,
         #POS_IF_FLU, 
         #TP_FLU_IF,
         PCR_RESUL,
         TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT, PCR_FLUBLI, FLUBLI_OUT,CLASSI_FIN,
         VACINA,
         EVOLUCAO, DT_EVOLUCA,VACINA_COV)

flu2023 <- 
  in2023_2 %>% 
  mutate(date_bth = dmy(DT_NASC),
         date_flu = dmy(DT_SIN_PRI),
         date_not = dmy(DT_NOTIFIC),
         date_dth = dmy(DT_EVOLUCA),
         sex = CS_SEXO %>% str_to_lower(),
         year = year(date_flu),
         
         typ = case_when(TP_FLU_PCR == 1 ~ "a",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         sub = case_when(PCR_FLUASU == 1 ~ "h1",
                         PCR_FLUASU == 2 ~ "h3",
                         PCR_FLUASU %in% 3:6 ~ "an",
                         FLUASU_OUT == 1 ~ "ao",
                         TP_FLU_PCR == 2 ~ "b",
                         TRUE ~ "other"),
         antigen=case_when(TP_FLU_AN == 1 ~"a",
                           TP_FLU_AN == 2 ~"b",
                           TRUE ~ "other"),
         classi_fin = case_when(CLASSI_FIN == 1 ~ "influenza",
                                CLASSI_FIN == 2 ~ "other_resp",
                                CLASSI_FIN == 3 ~ "other_etiol",
                                CLASSI_FIN == 4 ~ "unk",
                                CLASSI_FIN == 5 ~ "covid",
                                TRUE ~ "unk"),
         # revise type and sub to "ab" when classi_fin==influenza
         typ = case_when((typ == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ typ),
         sub = case_when((sub == "other" & classi_fin == "influenza") ~ "ab",
                         TRUE ~ sub),
         flu = ifelse(POS_PCRFLU == 1 | 
                        typ != "other" | 
                        sub != "other", 
                      1, 0),
         typ = case_when((typ=="other" & flu==0 & classi_fin!="covid") ~ "notflu",
                         TRUE ~typ),
         sub = case_when((sub=="other" & flu==0 & classi_fin!="covid") ~ "notflu",
                         TRUE ~sub),
         outcome = case_when(EVOLUCAO == 1 ~ "survived",
                             EVOLUCAO == 2 ~ "death_flu",
                             EVOLUCAO == 3 ~ "death_oth",
                             EVOLUCAO == 4 ~ "death_inv",
                             EVOLUCAO == 5 ~ "missing",
                             TRUE ~ "oth"),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk"),
         pcr_resul = case_when(PCR_RESUL == 1 ~ "Detectable", 
                               PCR_RESUL == 2 ~ "Not Detectable", 
                               PCR_RESUL == 3 ~ "Inconclusive", 
                               PCR_RESUL == 4 ~ "Not Realized", 
                               PCR_RESUL == 5 ~ "Waiting for Result", 
                               PCR_RESUL == 9 ~ "Ignored",
                               TRUE ~ "unk"),
         vaccine_cov = case_when(VACINA_COV == 1 ~ "y",
                                VACINA_COV == 2 ~ "n",
                                VACINA_COV == 9 ~ "unk",
                                TRUE ~ "na"))%>% 
  #filter(flu == 1) %>% 
  select(year, sex, date_bth, date_flu, date_not, date_dth, 
         typ, sub, classi_fin,
         outcome,vaccine)


flu2023 %>% 
  group_by(sub) %>% 
  summarise(n = n())

flu2023 %>% 
  group_by(outcome) %>% 
  summarise(n = n())


flu2023 %>% 
  group_by(vaccine) %>% 
  summarise(n = n())



# ~~~~~~~~~~~~~~~~~~~~~~
# All data together ====
# ~~~~~~~~~~~~~~~~~~~~~~
flu <- 
  bind_rows(flu0911, flu1218, flu19, flu2023)

write_rds(flu, "data_inter/flu_data_brazil_2009_2023.rds")

flu %>% 
  group_by(classi_fin) %>% 
  summarise(n = n())

flu %>% 
  group_by(sub) %>% 
  summarise(n = n())

