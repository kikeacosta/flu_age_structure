source("code/00_functions.R")

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

cut_db1 <- function(db){
  db2 <- 
    db %>% 
    mutate(id = 1:n()) %>% 
    select(
      id,
      DT_NOTIFIC,
      PCR_RES,
      PCR_ETIOL, PCR_TIPO_H, PCR_TIPO_N,
      CULT_RES,
      HEMA_RES, HEMA_ETIOL, HEM_TIPO_H, HEM_TIPO_N,
      CLASSI_FIN,
      PCR,
      RES_FLUA, RES_FLUB, RES_FLUASU, 
      # DS_OUTSUB,
      DT_SIN_PRI,
      DT_NOTIFIC,
      DT_NASC,
      CS_SEXO,
      NU_IDADE_N,
      HOSPITAL,
      EVOLUCAO,
    )
}
cut_db2 <- function(db){
  db2 <- 
    db %>% 
    mutate(id = 1:n()) %>% 
    select(
      id,
      DT_NOTIFIC,
      POS_IF_FLU, TP_FLU_IF,
      PCR_RESUL, POS_PCRFLU, 
      TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT,
      PCR_FLUBLI, FLUBLI_OUT 
    )
}

dt <-
  bind_rows(
    cut_db1(in09),
    cut_db1(in10),
    cut_db1(in11),
    cut_db1(in12),
    cut_db1(in13),
    cut_db1(in14),
    cut_db1(in15),
    cut_db1(in16),
    cut_db1(in17),
    cut_db1(in18)
  ) %>% 
  mutate(across(3:16, ~replace_na(.x, 0))) 

dt2 <- 
  dt %>% 
  mutate(
    # identifying influenza infections
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pcr = ifelse(PCR_RES %in% 1:3 | PCR == 1, 1, 0),
    pcr = ifelse(is.na(PCR_RES) & is.na(PCR), 0, pcr),
    test1 = ifelse(PCR_ETIOL == 1 | 
                     PCR_TIPO_H %in% 1:16 |
                     PCR_TIPO_N %in% 1:9, 
                   1, 0),
    test2 = ifelse(HEMA_RES == 1 | 
                     HEMA_ETIOL %in% 1:4 | 
                     HEM_TIPO_H %in% 1:16 |
                     HEM_TIPO_N %in% 1:9, 
                   1, 0),
    test3 = ifelse(CULT_RES == 1, 1, 0),
    test4 = ifelse(RES_FLUA == 1 | 
                     RES_FLUB == 1,
                   1, 0),
    test5 = ifelse(CLASSI_FIN == 1, 
                   1, 0),
    test = ifelse(test1 + test2 + test3 + test4 + test5 > 0,
                  1, 0),
    # identifying virus type (A/B) 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    typ1 = case_when(PCR_ETIOL %in% c(1, 2, 4) ~ "a",
                     PCR_ETIOL == 3 ~ "b",
                     PCR_TIPO_H %in% 1:6 ~ "a",
                     PCR_TIPO_N %in% 1:9 ~ "a",
                     TRUE ~ "z"),
    typ2 = case_when(HEMA_ETIOL %in% c(1, 2, 4) ~ "a",
                     HEMA_ETIOL %in% c(3) ~ "b",
                     HEM_TIPO_H %in% 1:16 ~ "a",
                     HEM_TIPO_N %in% 1:9 ~ "a",
                     TRUE ~ "z"),
    typ3 = case_when(RES_FLUA == 1 ~ "a",
                     RES_FLUASU %in% 1:6 ~ "a",
                     RES_FLUB == 1 ~ "b",
                     TRUE ~ "z"),
    typ = case_when(typ1 == "a" | typ2 == "a" | typ3 == "a" ~ "a",
                    typ1 == "b" | typ2 == "b" | typ3 == "b" ~ "b",
                    TRUE ~ "z"),
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # identifying A virus subtype (h1/h3)
    sub1 = case_when(PCR_TIPO_H == 1 ~ "h1",
                     PCR_TIPO_H == 3 ~ "h3",
                     is.na(PCR_TIPO_H) ~ "z",
                     TRUE ~ "z"),
    sub2 = case_when(HEM_TIPO_H == 1 ~ "h1",
                     HEM_TIPO_H == 3 ~ "h3",
                     TRUE ~ "z"),
    sub3 = case_when(RES_FLUASU == 1 ~ "h1",
                     RES_FLUASU == 2 ~ "h1",
                     RES_FLUASU == 3 ~ "h3",
                     RES_FLUASU == 5 ~ "h3",
                     TRUE ~ "z"),
    sub = case_when(sub1 == "h1" | sub2 == "h1" | sub3 == "h1" ~ "h1",
                    sub1 == "h3" | sub2 == "h3" | sub3 == "h3" ~ "h3",
                    TRUE ~ "z")
  ) %>% 
  mutate(flu = ifelse(test == 1 | 
                        typ %in% c("a", "b") | 
                        sub %in% c("h1", "h3"), 
                      1, 0)) %>% 
  mutate(
    outcome = case_when(EVOLUCAO == 1 ~ "survived",
                        EVOLUCAO == 2 ~ "death_flu",
                        EVOLUCAO == 3 ~ "death_oth",
                        EVOLUCAO == 4 ~ "death_inv",
                        EVOLUCAO == 5 ~ "missing",
                        TRUE ~ "oth"),
    hosp = ifelse(HOSPITAL == 1, 1, 0),
    date1 = dmy(DT_SIN_PRI),
    date2 = dmy(DT_NOTIFIC),
    date_bth = dmy(DT_NASC),
    year = ifelse(is.na(date1),
                  year(date2), 
                  year(date1)),
    sex = CS_SEXO %>% str_to_lower(),
    age = case_when(NU_IDADE_N < 4000 ~ 0,
                    NU_IDADE_N > 4000 ~ NU_IDADE_N - 4000),
    cohort = ifelse(is.na(date_bth),
                    year - age, 
                    year(date_bth)),
  )


flu0918 <- 
  dt2 %>%  
  select(year, sex, age, cohort, flu, typ, sub, hosp, outcome) %>% 
  drop_na()

table(flu0918$year)
table(flu0918$year, flu0918$flu)
table(flu0918$year, flu0918$typ)
table(flu0918$year, flu0918$sub)
table(flu0918$year, flu0918$outcome)

all <- 
  flu0918 %>% 
  summarise(n = n(), .by = cohort) %>% 
  mutate(cx = n / sum(n), 
         status = "all SRAG")

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
         hosp = ifelse(HOSPITAL == 1, 1, 0),
         date_dth = dmy(DT_EVOLUCA),
         sex = CS_SEXO %>% str_to_lower(),
         age = NU_IDADE_N,
         year = year(date_flu),
         cohort = ifelse(is.na(date_bth),
                         year - age, 
                         year(date_bth)),
         vaccine = case_when(VACINA == 1 ~ "y",
                             VACINA == 2 ~ "n",
                             VACINA == 9 ~ "unk",
                             TRUE ~ "unk")
         )


flu19 <- 
  in193 %>% 
  # filter(flu == 1) %>% 
  select(year, sex, age, cohort, flu, typ, sub, hosp, outcome)

# ~~~~~~~~~~~~~~~~~~~~~~
# All data together ====
# ~~~~~~~~~~~~~~~~~~~~~~
flu <- 
  bind_rows(flu0918, flu19)

test <- 
  flu %>% 
  group_by(year, age) %>% 
  summarise(n = n())

test %>% 
  spread(year, n)


write_rds(flu, "data_inter/flu_data_brazil_2009_2019_v3.rds",
          compress = "xz")

# comparison with previous versions
flu_prev1 <- read_rds("data_inter/flu_data_brazil_2009_2019.rds")
flu_prev2 <- read_rds("data_inter/flu_data_brazil_2009_2019_v2.rds")

flu_prev1 %>% 
  group_by(sub) %>% 
  summarise(n = n()) %>%  
  ungroup()
flu_prev2 %>% 
  group_by(sub) %>% 
  summarise(n = n()) %>%  
  ungroup()
flu %>% 
  group_by(sub) %>% 
  summarise(n = n()) %>%  
  ungroup()

# write_rds(flu_prev1, "data_inter/flu_data_brazil_2009_2019.rds",
#           compress = "xz")
# write_rds(flu_prev2, "data_inter/flu_data_brazil_2009_2019_v2.rds",
#           compress = "xz")


















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
