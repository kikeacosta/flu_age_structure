rm(list=ls())
source("code/00_functions.R")

epi <- read_csv("data_input/VIW_FID_EPI.csv")

epi2 <- 
  epi %>% 
  select(code = COUNTRY_CODE, 
         date = ISO_WEEKSTARTDATE,
         year = ISO_YEAR,
         week = ISO_WEEK,
         age_label = AGEGROUP_CODE,
         ili_css = ILI_CASES,
         ili_pts = ILI_OUTPATIENTS,
         ili_pop_cov = ILI_POP_COV, 
         sari_css = SARI_CASES, 
         sari_pts = SARI_INPATIENTS)


flu <- read_csv("data_input/VIW_FNT.csv")

flu2 <- 
  flu %>% 
  select(code = COUNTRY_CODE, 
         source = ORIGIN_SOURCE,
         date = ISO_WEEKSTARTDATE,
         year = ISO_YEAR,
         week = ISO_WEEK,
         h1 = AH1N12009,
         h1pre = AH1, 
         h3 = AH3,
         inf_a = INF_A,
         inf_b = INF_B,
         inf = INF_ALL,
         inf_neg = INF_NEGATIVE,
         ade = ADENO,
         boc = BOCA,
         cor = HUMAN_CORONA,
         met = METAPNEUMO,
         par = PARAINFLUENZA,
         rhn = RHINO,
         rsv = RSV,
         oth = OTHERRESPVIRUS,
         ili = ILI_ACTIVITY) %>% 
  gather(-code, -source, -date, -year, -week, key = measure, value = value) %>% 
  replace_na(list(value = 0)) %>% 
  spread(measure, value) %>% 
  mutate(
    # # adjusting for flu infections
    inf = ifelse(inf > inf_a + inf_b, inf, inf_a + inf_b), 
    # number of samples with positive
    vir = inf + rsv + ade + boc + cor + met + par + rhn + oth,
    # number of H3 and H1 infections
    h1h3 = h1 + h3,
    p_h1 = h1 / vir,
    p_h3 = h3 / vir,
    p_inf = inf / vir)

flu3 <- 
  flu2 %>% 
  filter(p_h1 >= 0.7 | p_h3 >= 0.7,
         h1h3 >= 200) %>% 
  select(code, year, week, inf, h1, h3, p_h1, p_h3) %>% 
  unique()

sel_combs <- 
  flu3 %>% 
  select(code, year, week) %>% 
  unique()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ilis <- 
  epi2 %>% 
  inner_join(sel_combs)

unique(ilis$code)

ilis2 <- 
  ilis %>% 
  mutate(age_lab2 = age_label) %>% 
  separate(age_lab2, c("age", "age_up"), sep = "TO") %>% 
  mutate(age = age %>% as.double(),
         age_up = age_up %>% as.double(),
         age = ifelse(age < 1, age * 10, age),
         age_up = ifelse(age_up < 1, age_up * 10, age_up),
         age_spn = age_up - age + 1) %>% 
  replace_na(list(ili_css = 0,
                  ili_pts = 0,
                  ili_pop_cov = 0, 
                  sari_css = 0, 
                  sari_pts = 0)) %>% arrange(code, date) %>% 
  group_by(code, date) %>% 
  mutate(age_grs = n()) %>% 
  ungroup() %>% 
  mutate(ili = ili_css + ili_pts,
         sari = sari_css + sari_pts) %>% 
  arrange(code, date, age)


