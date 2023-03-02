rm(list=ls())
source("code/00_functions.R")


epi <- read_csv("data_input/VIW_FID_EPI.csv")

epi2 <- 
  epi %>% 
  select(code = COUNTRY_CODE, 
         date = ISO_WEEKSTARTDATE,
         year = ISO_YEAR,
         week = ISO_WEEK,
         age = AGEGROUP_CODE,
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
         h5 = AH5,
         h7 = AH7N9,
         a_non_sub = ANOTSUBTYPED,
         a_imp_sub = ANOTSUBTYPABLE,
         a_other = AOTHER_SUBTYPE,
         inf_a = INF_A,
         inf_b = INF_B,
         rsv = RSV,
         inf_all = INF_ALL,
         css_neg = INF_NEGATIVE,
         css = SPEC_PROCESSED_NB,
         ili = ILI_ACTIVITY) %>% 
  gather(-code, -source, -date, -year, -week, key = measure, value = value) %>% 
  replace_na(list(value = 0)) %>% 
  spread(measure, value)

flu3 <- 
  flu2 %>% 
  mutate(h1h3 = h1+h3,
         css_pos = css - css_neg,
         ph1 = ifelse(css_pos > 0, h1 / css_pos, 0),
         ph3 = ifelse(css_pos > 0, h3 / css_pos, 0),
         ph1h3 = h1h3 / css_pos)

# %>% 
#   mutate(css_all = css_neg + inf_all)

# test <- 
#   flu_sub %>% 
#   group_by(code, source, date, year, week, measure) %>% 
#   summarise(n = n()) %>% 
#   ungroup()



# %>% 
#   spread(measure, value)
#   mutate(
#     ph10 = h1 / (h1 + h3),
#     ph12 = h1 / inf_a,
#     ph13 = h1 / all_inf,
#     ph33 = h3 / all_inf)

test1 <- 
  flu_sub %>% 
  filter(all_inf != inf_a + inf_b)

test <- 
  flu_sub %>% 
  filter((ph13 > 0.9 & ph13 <= 1) | (ph33 > 0.9 & ph33 <= 1))


# 





     