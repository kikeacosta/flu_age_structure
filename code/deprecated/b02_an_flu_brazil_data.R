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
  summarise(css = n()) %>% 
  ungroup() 

# 
flu0 %>% 
  filter(sub %in% c("h1", "h3")) %>% 
  ggplot()+
  geom_line(aes(date_flu, css, col = sub))+
  theme_bw()


# flu circulation by sex, age, and subtype ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flu2 <- 
  flu %>% 
  mutate(cohort = year(date_bth),
         age = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
         age = round(age)) %>% 
  group_by(year, sex,  age, typ, sub) %>% 
  summarise(css = n()) %>% 
  ungroup()

flu3 <- 
  flu2 %>% 
  group_by(year, age, typ, sub) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(flu2) %>% 
  mutate(cohort = year - age) %>% 
  left_join(pop3) %>% 
  drop_na(age) %>% 
  mutate(ins = 1e5*css / pop)

flu3 %>% 
  group_by(sub) %>% 
  summarise(sum(css))

write_rds(flu3, "data_inter/flu_data_brazil_exposures_2009_2019.rds")

