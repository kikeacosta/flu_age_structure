rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
library(wpp2022)

options(scipen=999)

flu <- 
  read_rds("data_inter/flu_data_brazil_2009_2019_v3.rds")

unique(flu$year)

# exposures
data(popAge1dt)

pop <- 
  popAge1dt %>% 
  select(year, name, age, m = popM, f = popF, t = pop) %>% 
  filter(name == "Brazil",
         year %in% 2000:2022) %>% 
  pivot_longer(c(m, f, t), 
               names_to = "sex", values_to = "pop") %>% 
  select(year, age, sex, pop) %>% 
  mutate(pop = pop*1e3)

unique(pop$sex)

# pop_f <- read_delim("data_input/brazil/pop_female.csv",
#                     skip = 4,
#                     delim = ";")
# 
# pop_m <- read_delim("data_input/brazil/pop_male.csv",
#                     skip = 4,
#                     delim = ";")
# 
# pop <- 
#   bind_rows(pop_f %>% mutate(sex = "f"),
#             pop_m %>% mutate(sex = "m")) %>% 
#   rename(age = 1) %>% 
#   gather(-sex, -age, key = year, value = pop) %>% 
#   drop_na(pop) %>% 
#   mutate(year = year %>% as.integer()) %>% 
#   separate(age, c("age", "trash")) %>% 
#   filter(age != "Total") %>% 
#   mutate(age = age %>% as.integer()) %>% 
#   select(-trash)
# 
# pop2 <- 
#   pop %>% 
#   reframe(pop = sum(pop), .by = c(age, year)) %>% 
#   mutate(sex = "t") %>% 
#   bind_rows(pop) %>% 
#   select(year, sex, age, pop) %>% 
#   arrange(year, sex, age)
# 
# unique(pop2$year)

# ungrouping exposures at ages 90-100
chunk <- 
  pop %>% 
  filter(sex == "t",
         year == 2010)

ungroup_open_age <- 
  function(chunk){
    
    ages <- chunk$age
    pops <- chunk$pop
    
    nlast <- 111 - max(ages)
    yr <- unique(chunk$year)
    sx <- unique(chunk$sex)
    
    V1 <- pclm(x = ages,
               y = pops,
               nlast = nlast)$fitted
    
    new_pop <- tibble(year = yr,
                      sex = sx,
                      age = 0:110, 
                      pop = V1)
  
    chunk %>% 
      filter(age < 90) %>% 
      bind_rows(new_pop %>% filter(age >= 90))
  }

pop2 <- 
  pop %>% 
  group_by(year, sex) %>% 
  do(ungroup_open_age(chunk = .data)) %>% 
  ungroup() 

# total flu circulation ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
flu0 <- 
  flu %>% 
  summarise(value = n(), .by = c(cohort, sub))

# 
flu0 %>% 
  filter(sub %in% c("h1", "h3")) %>% 
  ggplot()+
  geom_line(aes(cohort, value, col = sub))+
  theme_bw()

# flu circulation by sex, age, and subtype ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flu2 <- 
  flu %>% 
  # filter(sub %in% c("h1", "h3", "b")) %>% 
  # mutate(age2 = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
  #        age2 = round(age2),
  #        age = ifelse(is.na(age), age2, age)) %>% 
  # select(-age2) %>% 
  summarise(value = n(), .by = c(year, sex, age, sub, hosp, outcome)) %>% 
  drop_na(age)

css <- 
  flu2 %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "cases")

hsp <- 
  flu2 %>% 
  filter(hosp == 1) %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "hosps")

dts <- 
  flu2 %>% 
  filter(outcome == "death_flu") %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "deaths")

flu3 <- 
  bind_rows(dts, css) %>% 
  group_by(year, age, outcome, sub) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(dts, css) %>% 
  mutate(cohort = year - age) %>% 
  left_join(pop2) %>% 
  drop_na(age) %>% 
  mutate(mx = 1e5 * value / pop)

flu3 %>% 
  filter(sex == "t") %>% 
  group_by(sub, outcome) %>% 
  summarise(sum(value))

write_rds(flu3, "data_inter/flu_data_brazil_exposures_2009_2019.rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flu3 <- read_rds("data_inter/flu_data_brazil_exposures_2009_2019.rds")

# H1 ====
# ~~~~~~~
h1 <- 
  flu3 %>% 
  filter(sub == "h1",
         sex == "t",
         outcome == "deaths") %>% 
  group_by(cohort, sex) %>% 
  summarise(css = sum(css),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(ins = 1e5 * css / pop)

h1 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, css))+
  geom_line(aes(cohort, css))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()+
  labs(title = "H1")

h1 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, ins))+
  geom_line(aes(cohort, ins))+
  geom_vline(xintercept = c(1918, 1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "H1")


# H3 ====
# ~~~~~~~
h3 <- 
  flu3 %>% 
  filter(sub == "h3",
         sex == "t",
         outcome == "death") %>% 
  group_by(cohort, sex) %>% 
  summarise(css = sum(css),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(ins = 1e5 * css / pop)

h3 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, css))+
  geom_line(aes(cohort, css))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()+
  labs(title = "H3")

h3 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, ins))+
  geom_line(aes(cohort, ins))+
  geom_vline(xintercept = c(1918, 1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "H3")

hs <- 
  flu3 %>% 
  filter(sub %in% c("h1", "h3"), sex == "t") %>% 
  group_by(sub, cohort) %>% 
  summarise(css = sum(css),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(mx = 1e5*css/pop) %>% 
  group_by(sub) %>% 
  mutate(cx_cs = css / sum(css),
         cx_mx = mx / sum(mx)) %>% 
  ungroup()

hs %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, cx_cs, col = sub))+
  geom_line(aes(cohort, cx_cs, col = sub))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()+
  labs(title = "cohort structure of infections")

hs %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, cx_mx, col = sub))+
  geom_line(aes(cohort, cx_mx, col = sub))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "cohort structure of infection rates")

