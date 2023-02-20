rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

# loading data ====
# ~~~~~~~~~~~~~~~~~.


# weekly positive influenza cases by subtype, sex, and age
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cs <- read_rds("data_inter/weekly_flu_subtypes_australia.rds")


# population
pop <- read_csv("data_input/Australia_exposure_hmd.csv")


# preparing population ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# only total sex, since 2000 and closing at age 100+

pop2 <- 
  pop %>% 
  filter(Year >= 2000) %>% 
  gather(-Year, -Age, key = sex, value = pop) %>% 
  rename(year = Year, age = Age) %>% 
  mutate(sex = sex %>% str_sub(1, 1) %>% str_to_lower(),
         age = ifelse(age == "110+", "110", age),
         age = age %>% as.integer(),
         age = age - age%%5,
         age = ifelse(age > 85, 85, age)) %>% 
  group_by(year, age, sex) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()


# annual configuration
cs2 <- 
  cs %>% 
  mutate(year = year(date),
         age = str_sub(age, 1, 2),
         age = age %>% as.integer()) %>% 
  group_by(year, sex, age, sub) %>% 
  summarise(css = sum(css)) %>% 
  ungroup()

# merging with population
# ~~~~~~~~~~~~~~~~~~~~~~~

cs3 <- 
  cs2 %>% 
  left_join(pop2) %>% 
  mutate(css = round(css))

write_rds(cs3, "data_inter/annual_flu_subtypes_australia.rds")



