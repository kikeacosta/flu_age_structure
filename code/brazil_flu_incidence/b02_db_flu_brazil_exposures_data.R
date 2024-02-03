rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

flu <- 
  read_rds("data_inter/flu_data_brazil_2009_2019.rds")

unique(flu$date_flu)

pop_f <- read_delim("data_input/brazil/pop_female.csv",
                    skip = 4,
                    delim = ";")

pop_m <- read_delim("data_input/brazil/pop_male.csv",
                    skip = 4,
                    delim = ";")

pop <- 
  bind_rows(pop_f %>% mutate(sex = "f"),
            pop_m %>% mutate(sex = "m")) %>% 
  rename(age = 1) %>% 
  gather(-sex, -age, key = year, value = pop) %>% 
  drop_na(pop) %>% 
  mutate(year = year %>% as.integer()) %>% 
  separate(age, c("age", "trash")) %>% 
  filter(age != "Total") %>% 
  mutate(age = age %>% as.integer()) %>% 
  select(-trash)

pop2 <- 
  pop %>% 
  group_by(age, year) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(pop) %>% 
  select(year, sex, age, pop) %>% 
  arrange(year, sex, age)


# ungrouping exposures at ages 90-100
chunk <- 
  pop2 %>% 
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

pop3 <- 
  pop2 %>% 
  group_by(year, sex) %>% 
  do(ungroup_open_age(chunk = .data)) %>% 
  ungroup() 

# total flu circulation ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
flu0 <- 
  flu %>% 
  mutate(cohort = year(date_bth),
         age = year - cohort) %>% 
  group_by(date_flu, sub) %>% 
  summarise(value = n()) %>% 
  ungroup() 

# 
flu0 %>% 
  filter(sub %in% c("h1", "h3")) %>% 
  ggplot()+
  geom_line(aes(date_flu, value, col = sub))+
  theme_bw()


# flu circulation by sex, age, and subtype ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flu2 <- 
  flu %>% 
  filter(sub %in% c("h1", "h3", "b")) %>% 
  mutate(age2 = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
         age2 = round(age2),
         age = ifelse(is.na(age), age2, age)) %>% 
  select(-age2) %>% 
  group_by(year, sex, age, sub, outcome) %>% 
  summarise(value = n()) %>% 
  ungroup() %>% 
  drop_na(age)

css <- 
  flu2 %>% 
  group_by(year, sex, age, sub) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(outcome = "cases")

dts <- 
  flu2 %>% 
  filter(outcome == "death_flu") %>% 
  group_by(year, sex, age, sub) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(outcome = "deaths")


flu3 <- 
  bind_rows(dts, css) %>% 
  group_by(year, age, outcome, sub) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(dts, css) %>% 
  mutate(cohort = year - age) %>% 
  left_join(pop3) %>% 
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

