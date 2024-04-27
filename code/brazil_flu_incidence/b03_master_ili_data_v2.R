rm(list=ls())
source("code/brazil_flu_incidence/b00_functions.R")

# all Brazil included in the data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop <- 
  read_rds("data_inter/brazil_exposures_1999_2022_ages_0_100.rds")
  
unique(pop$age)
unique(pop$sex)

dt <- 
  read_rds("data_inter/flu_data_brazil_2009_2023_states.rds") %>% 
  mutate(age = year - cohort) %>% 
  filter(sex != 1) 


table(dt$typ, dt$sub)
unique(dt2$iso2)
# keeping only subtype and redistributing deaths by sex
dt2 <- 
  dt %>% 
  # filter (iso2 %in% c(states_sel)) %>% 
  mutate(age = year - cohort,
         age = case_when(age < 0 ~ 0, 
                         age > 100 ~ 100,
                         TRUE ~ age)) %>% 
  select(-typ) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age, flu, sub, hosp, noflu, outcome)) %>% 
  # imputing sex distribution
  spread(sex, value) %>% 
  replace_na(list(f = 0, m = 0, i = 0)) %>% 
  mutate(
    t = f+m+i,
    dts_sum = f+m,
    f = round(f*t/dts_sum),
    m = round(m*t/dts_sum)
  ) %>% 
  select(-dts_sum, -i) %>% 
  gather(t, f, m, key = sex, value = value) %>% 
  mutate(dth = ifelse(outcome %in% c("death_flu", "death_inv", "death_oth"), 
                      1, 0))


# grouping data 
# ili
ili <- 
  dt2 %>%
  filter(sex == "t") %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "ili")

# sari
sari <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(hosp == 1) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "sari")

# flu cases
flu <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1) %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "flu")

# H1 cases
h1 <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1,
         sub == "h1") %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "h1")

# H3 cases
h3 <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1,
         sub == "h3") %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "h3")

non_flu <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(noflu == 1) %>% 
  # mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "nonflu")

# putting all together
all <- 
  bind_rows(ili, sari, flu, h1, h3, non_flu) %>% 
  complete(year, sex, age, type, 
           fill = list(css = 0, hsp = 0, dts = 0)) %>% 
  arrange(type, year, sex, age) %>% 
  left_join(pop2) 

unique(all$type)

write_rds(all, "data_inter/master_brazil_flu_2009_2023.rds")

