rm(list=ls())
source("code/00_functions.R")

# Brazil
# ~~~~~~
# 
# # to look at the cohort composition of each flu status
# library(wpp2022)
# data(popAge1dt)
# br_pop <- 
#   popAge1dt %>% 
#   filter(name == "Brazil",
#          year %in% c(2010, 2016, 2019)) %>% 
#   mutate(cohort = year - age) %>% 
#   select(year, cohort, n = pop) %>% 
#   group_by(year) %>% 
#   mutate(cx = n/sum(n),
#          status = "all brazil wpp")
# 
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
in20 <- read_delim("data_input/brazil/INFLUD20.csv", delim = ";")
in21 <- read_delim("data_input/brazil/INFLUD21.csv", delim = ";")
in22 <- read_delim("data_input/brazil/INFLUD22.csv", delim = ";")
in23 <- read_delim("data_input/brazil/INFLUD23.csv", delim = ";")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unique(in09$NU_IDADE_N) %>% sort()
unique(in11$NU_IDADE_N) %>% sort()
unique(in18$NU_IDADE_N) %>% sort()
unique(in19$NU_IDADE_N) %>% sort()
unique(in21$NU_IDADE_N) %>% sort()

# states
st <- read_csv("data_input/brazil/states_codes.csv", 
               locale = readr::locale(encoding = "latin1"))

unique(in09$SG_UF_NOT) %>% sort

# functions to format data inputs ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2009-2018
cut_db1 <- 
  function(db){
  db %>% 
    mutate(id = 1:n()) %>% 
    select(
      id,
      SG_UF_NOT,
      DT_NOTIFIC,
      DT_SIN_PRI,
      DT_NOTIFIC,
      DT_NASC,
      CS_SEXO,
      NU_IDADE_N,
      PCR_RES,
      PCR_ETIOL, PCR_TIPO_H, PCR_TIPO_N,
      CULT_RES,
      HEMA_RES, HEMA_ETIOL, HEM_TIPO_H, HEM_TIPO_N,
      CLASSI_FIN,
      PCR,
      RES_FLUA, RES_FLUB, RES_FLUASU, 
      # DS_OUTSUB,
      HOSPITAL,
      EVOLUCAO,
    )
}

# 2019
cut_db2 <- 
  function(db){
  db %>% 
    mutate(id = 1:n()) %>% 
    select(
      id,
      SG_UF_NOT,
      DT_NOTIFIC,
      DT_NASC,
      CS_SEXO,
      NU_IDADE_N,
      DT_SIN_PRI,
      
      POS_IF_FLU, TP_FLU_IF,
      PCR_RESUL, POS_PCRFLU, 
      TP_FLU_PCR, PCR_FLUASU, FLUASU_OUT,
      PCR_FLUBLI, FLUBLI_OUT,
      CLASSI_FIN,
      HOSPITAL,
      EVOLUCAO,
    )
}

# 2020-2023
cut_db3 <- 
  function(db){
  db %>% 
    mutate(id = 1:n()) %>% 
    select(
      id,
      SG_UF_NOT,
      DT_NOTIFIC,
      DT_NASC,
      CS_SEXO,
      NU_IDADE_N,
      DT_SIN_PRI,
      
      # etiologic diagnosis of flu
      POS_AN_FLU,
      # A/B
      TP_FLU_AN,
      # PCR
      PCR_RESUL,
      # PCR flu
      POS_PCRFLU, 
      TP_FLU_PCR, 
      # A
      PCR_FLUASU, FLUASU_OUT,
      # B
      PCR_FLUBLI, FLUBLI_OUT,

      CLASSI_FIN,
      HOSPITAL,
      EVOLUCAO,
    )
  }



# ~~~~~~~~~~~~~~~~~~~~~
# data 2009 - 2018 ====
# ~~~~~~~~~~~~~~~~~~~~~

dt1 <-
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
  mutate(across(7:22, ~replace_na(.x, 0))) 

# table(dt1$PCR_RES)
# table(dt1$PCR_ETIOL)
# table(dt1$HEMA_RES)
# table(dt1$CULT_RES)
# table(dt1$SG_UF_NOT)

flu0918 <- 
  dt1 %>% 
  mutate(
    # identifying influenza infections
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pcr = ifelse(PCR_RES %in% 1:3 | PCR == 1, 1, 0),
    pcr = ifelse(is.na(PCR_RES) & is.na(PCR), 0, pcr),
    test1 = ifelse(PCR_ETIOL %in% 1:4 | 
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
                     RES_FLUASU %in% 1:6 | 
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
    typ = ifelse(typ == "z" & (CLASSI_FIN == 1 | test == 1), "ab", typ),
    # identifying A virus subtype (h1/h3)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  mutate(flu = ifelse(CLASSI_FIN == 1 |
                        test == 1 | 
                        typ %in% c("a", "b", "ab") | 
                        sub %in% c("h1", "h3"), 
                      1, 0),
         # h1 = ifelse(sub == "h1", 1, 0),
         # h3 = ifelse(sub == "h3", 1, 0),
         noflu = ifelse(RES_FLUA == 2 & RES_FLUB == 2, 1, 0)) %>% 
  mutate(
    outcome = case_when(EVOLUCAO == 1 ~ "survived",
                        EVOLUCAO == 2 ~ "death_flu",
                        EVOLUCAO == 3 ~ "death_oth",
                        EVOLUCAO == 4 ~ "death_inv",
                        EVOLUCAO == 5 ~ "missing",
                        TRUE ~ "missing"),
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
    age2 = interval(ymd(date_bth), ymd(date1)) %>% as.numeric('years'),
    age2 = round(age2),
    age = ifelse(is.na(age), age2, age),
    cohort = ifelse(is.na(date_bth),
                    year - age, 
                    year(date_bth)),
  ) %>% 
  select(cod = SG_UF_NOT, year, sex, age, cohort, 
         flu, typ, sub, hosp, noflu, outcome) %>%  
  left_join(st %>% select(-state)) %>% 
  drop_na() %>% 
  select(-cod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~
# data 2019 ====
# ~~~~~~~~~~~~~~

dt2 <-
  bind_rows(
    cut_db2(in19),
  ) %>% 
  mutate(across(c(8:13, 15, 17:19), ~replace_na(.x, 0)))


flu19 <- 
  dt2 %>% 
  mutate(
    test1 = case_when(POS_PCRFLU == 1 ~ "pos",
                      POS_PCRFLU == 2 ~ "neg",
                      POS_PCRFLU == 3 ~ "inc",
                      PCR_RESUL == 4 ~ "not",
                      TRUE ~ "z"),
    test2 = case_when(POS_IF_FLU == 1 ~ "pos",
                      POS_IF_FLU == 2 ~ "neg",
                      POS_IF_FLU == 3 ~ "inc",
                      TRUE ~ "z"),
    typ1 = case_when(TP_FLU_PCR == 1 ~ "a",
                     TP_FLU_PCR == 2 ~ "b",
                     TRUE ~ "z"),
    typ2 = case_when(TP_FLU_IF == 1 ~ "a",
                     TP_FLU_IF == 2 ~ "b",
                     TRUE ~ "z"),
    sub1 = case_when(PCR_FLUASU == 1 ~ "h1",
                     PCR_FLUASU == 2 ~ "h3",
                     PCR_FLUASU %in% 3:5 ~ "an",
                     PCR_FLUASU == 6 ~ "h3",
                     TP_FLU_PCR == 2 ~ "b",
                     TRUE ~ "z"),
    sub2 = case_when(str_detect(FLUASU_OUT, "H1") ~ "h1",
                     str_detect(FLUASU_OUT, "H3") ~ "h3",
                     TRUE ~ "z")
  ) %>% 
  mutate(
    flu = ifelse(CLASSI_FIN == 1 |
                   test1 == "pos" | test2 == "pos" | 
                   typ1 %in% c("a", "b") | typ2 %in% c("a", "b") | 
                   sub1 %in% c("h1", "h3") | sub2 %in% c("h1", "h3"), 
                 1, 0),
    noflu = ifelse(test1 == 2, 1, 0),
    typ = case_when(typ1 == "a" | typ2 == "a" ~ "a",
                    typ1 == "b" | typ2 == "b" ~ "b",
                    TRUE ~ "z"),
    typ = ifelse(typ == "z" & flu == 1, "ab", typ),
    sub = case_when(sub1 == "h1" | sub2 == "h1" ~ "h1",
                    sub1 == "h3" | sub2 == "h3" ~ "h3",
                    TRUE ~ "z"),
    outcome = case_when(EVOLUCAO == 1 ~ "survived",
                        EVOLUCAO == 2 ~ "death_flu",
                        EVOLUCAO == 3 ~ "death_oth",
                        EVOLUCAO == 4 ~ "death_inv",
                        EVOLUCAO == 5 ~ "missing",
                        TRUE ~ "oth"),
    hosp = ifelse(HOSPITAL == 1, 1, 0)
  ) %>% 
  mutate(
    date_bth = dmy(DT_NASC),
    date_flu = dmy(DT_SIN_PRI),
    date_not = dmy(DT_NOTIFIC),
    # date_vac = dmy(DT_UT_DOSE),
    year = year(date_flu),
    sex = CS_SEXO %>% str_to_lower(),
    age = NU_IDADE_N,
    age2 = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
    age2 = round(age2),
    age = ifelse(is.na(age), age2, age),
    cohort = ifelse(is.na(date_bth),
                    year - age, 
                    year(date_bth))
  ) %>% 
  select(iso2 = SG_UF_NOT, year, sex, age, cohort, 
         flu, typ, sub, hosp, outcome)

table(flu19$sex)
table(flu0918$sex)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~
# data 2020 - 2023 ====
# ~~~~~~~~~~~~~~~~~~~~~

dt3 <-
  bind_rows(
    # cut_db2(in19),
    cut_db3(in20),
    cut_db3(in21),
    cut_db3(in22),
    cut_db3(in23)
  ) %>% 
  mutate(across(8:19, ~replace_na(.x, 0))) 


flu2023 <- 
  dt3 %>% 
  mutate(
    test1 = case_when(POS_PCRFLU == 1 ~ "pos",
                      POS_PCRFLU == 2 ~ "neg",
                      POS_PCRFLU == 3 ~ "inc",
                      PCR_RESUL == 4 ~ "not",
                      TRUE ~ "z"),
    test2 = case_when(POS_AN_FLU == 1 ~ "pos",
                      POS_AN_FLU == 2 ~ "neg",
                      POS_AN_FLU == 3 ~ "inc",
                      TRUE ~ "z"),
    typ1 = case_when(TP_FLU_PCR == 1 ~ "a",
                     TP_FLU_PCR == 2 ~ "b",
                     TRUE ~ "z"),
    typ2 = case_when(TP_FLU_AN == 1 ~ "a",
                     TP_FLU_AN == 2 ~ "b",
                     TRUE ~ "z"),
    sub1 = case_when(PCR_FLUASU == 1 ~ "h1",
                     PCR_FLUASU == 2 ~ "h3",
                     PCR_FLUASU %in% 3:5 ~ "an",
                     PCR_FLUASU == 6 ~ "h3",
                     TP_FLU_PCR == 2 ~ "b",
                     TRUE ~ "z"),
    sub2 = case_when(str_detect(FLUASU_OUT, "H1") ~ "h1",
                     str_detect(FLUASU_OUT, "H3") ~ "h3",
                     TRUE ~ "z")
  ) %>% 
  mutate(
    flu = ifelse(CLASSI_FIN == 1 |
                   test1 == "pos" | test2 == "pos" | 
                   typ1 %in% c("a", "b") | typ2 %in% c("a", "b") | 
                   sub1 %in% c("h1", "h3") | sub2 %in% c("h1", "h3"), 
                 1, 0),
    noflu = ifelse(test1 == 2, 1, 0),
    typ = case_when(typ1 == "a" | typ2 == "a" ~ "a",
                    typ1 == "b" | typ2 == "b" ~ "b",
                    TRUE ~ "z"),
    typ = ifelse(typ == "z" & flu == 1, "ab", typ),
    sub = case_when(sub1 == "h1" | sub2 == "h1" ~ "h1",
                    sub1 == "h3" | sub2 == "h3" ~ "h3",
                    TRUE ~ "z"),
    outcome = case_when(EVOLUCAO == 1 ~ "survived",
                        EVOLUCAO == 2 ~ "death_flu",
                        EVOLUCAO == 3 ~ "death_oth",
                        EVOLUCAO == 4 ~ "death_inv",
                        EVOLUCAO == 5 ~ "missing",
                        TRUE ~ "oth"),
    hosp = ifelse(HOSPITAL == 1, 1, 0)
  ) %>% 
  mutate(
    date_bth = dmy(DT_NASC),
    date_flu = dmy(DT_SIN_PRI),
    date_not = dmy(DT_NOTIFIC),
    # date_vac = dmy(DT_UT_DOSE),
    year = year(date_flu),
    sex = CS_SEXO %>% str_to_lower(),
    age = NU_IDADE_N,
    age2 = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
    age2 = round(age2),
    age = ifelse(is.na(age), age2, age),
    cohort = ifelse(is.na(date_bth),
                    year - age, 
                    year(date_bth))
  ) %>% 
  select(iso2 = SG_UF_NOT, year, sex, age, cohort, 
         flu, typ, sub, hosp, outcome)

table(flu19$sex)
table(flu0918$sex)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~
# All data together ====
# ~~~~~~~~~~~~~~~~~~~~~~
flu <- 
  bind_rows(flu0918, flu19, flu2023) %>% 
  # adjusting types and subtypes to identify better flu cases
  mutate(sub = case_when(typ == "b" ~ "b",
                         typ == "a" & sub %in% c("an", "ao", "z") ~ "hx",
                         flu == 1 & sub == "z" ~ "ab",
                         TRUE ~ sub),
         typ = ifelse(flu == 1 & typ == "z", "ab", typ))

table(flu$typ, flu$sub)
table(flu$flu, flu$sub)
table(flu$flu, flu$typ)

write_rds(flu, "data_inter/flu_data_brazil_2009_2023_states.rds",
          compress = "xz")

unique(flu$typ)


# analysis by state ====
# ~~~~~~~~~~~~~~~~~~~~~~
# cases by state
dt_state <- 
  dt2 %>% 
  mutate(dth = ifelse(outcome %in% c("death_flu", "death_oth", "death_inv"), 1, 0),
         dth_flu = ifelse(outcome %in% c("death_flu", "death_oth", "death_inv"), 1, 0)) %>% 
  group_by(SG_UF_NOT, year) %>% 
  summarise(css = n(),
            hsp = sum(hosp),
            flu = sum(flu),
            h1 = sum(h1),
            h3 = sum(h3),
            dts = sum(dth))



stt_codes <- 
  read_csv("data_input/brazil/states_codes.csv", 
           locale = readr::locale(encoding = "latin1")) %>% 
  select(cod, iso2)

pop_st <- 
  read.delim("data_input/brazil/bra_population_state_2010_2025.csv", 
             skip = 3, 
             encoding = "UTF-8", sep = ";") %>% 
  drop_na() %>% 
  rename(year = Ano) %>% 
  mutate(year = year %>% as.numeric()) %>% 
  gather(-year, key = iso2, value = pop) %>% 
  left_join(stt_codes) %>% 
  filter(year %in% 2010:2023) %>% 
  group_by(year, iso2, cod) %>% 
  summarise(pop = mean(pop))

pop_st2 <- 
  pop_st %>% 
  filter(year == 2010) %>% 
  mutate(year = 2009) %>% 
  bind_rows(pop_st)


dt_state2 <- 
  dt_state %>% 
  rename(cod = SG_UF_NOT) %>% 
  left_join(pop_st2) %>% 
  mutate(inc_sari = css / pop,
         inc_flu = flu / pop,
         cdr = 1e5 * dts / pop) %>% 
  group_by(iso2) %>% 
  mutate(cdr_av = mean(cdr))

dt_state2 %>% 
  mutate(state = paste0(cod, "_", iso2)) %>% 
  ggplot()+
  geom_point(aes(cdr, state, col = year))+
  geom_point(aes(cdr_av, state), col = "red")+
  theme_bw()

# ggsave("figures/cdr_by_state.png",
#        w = 6, h = 6)



