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
cs <- read_rds("data_inter/annual_flu_subtypes_australia.rds")

# population
pop <- read_csv("data_input/Australia_exposure_hmd.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# preparing population ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# population in single-year of age
pop2 <- 
  pop %>% 
  filter(Year >= 2008) %>% 
  gather(-Year, -Age, key = sex, value = pop) %>% 
  rename(year = Year, age = Age) %>% 
  mutate(sex = sex %>% str_sub(1, 1) %>% str_to_lower(),
         age = ifelse(age == "110+", "110", age),
         age = age %>% as.integer(),
         age = ifelse(age > 100, 100, age)) %>% 
  group_by(year, age, sex) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()

# adding 2021 with the 2020 population
pop3 <- 
  pop2 %>% 
  bind_rows(pop2 %>% 
              filter(year == 2020) %>% 
              mutate(year = 2021))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ungrouping data ====
# ~~~~~~~~~~~~~~~~~~~~

chunk <- 
  cs %>% 
  filter(year == 2011,
         sex == "t",
         sub == "h1")

# ungrouping ages in single-year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ungroup_this <- 
  function(chunk,
           lambda = 1e10){

    ages <- chunk$age
    # increasing numbers to avoid inconsistent outcomes (decreasing afterwards)
    cass <- (chunk$css + 0.001) * 1e5
    # closing at age 100+
    nlast <- 101 - max(ages)
    yr <- unique(chunk$year)
    sx <- unique(chunk$sex)
    pops <- 
      pop3 %>% 
      filter(year == yr,
             sex == sx) %>% 
      pull(pop)
    
    sb <- unique(chunk$sub)
      
    cat(paste0(yr, "_", sx, "_", sb, "\n"))
    
    V1 <- pclm(x = ages,
               y = cass,
               nlast = nlast)$fitted
    V2 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               control = list(lambda = lambda, deg = 3))$fitted
    V3 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               offset = pops)$fitted * pops
    V4 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               offset = pops,
               control = list(lambda = lambda, deg = 3))$fitted * pops
    
    out <- 
      tibble(age = 0:100, 
             cs1 = V1/1e5,
             cs2 = V2/1e5,
             cs3 = V3/1e5,
             cs4 = V4/1e5)
    return(out)
  }

unique(cs$sub)

unique(cs$sub)

cs2 <- 
  cs %>% 
  filter(sub %in% c("h1", "h3"),
         year >= 2009) %>% 
  arrange(year, sex, sub, age) %>% 
  group_by(year, sex, sub) %>% 
  do(ungroup_this(chunk = .data, lambda = 1e10)) %>% 
  ungroup() %>% 
  gather(cs1, cs2, cs3, cs4, key = method, value = css) 

cs2 %>% 
  filter(sub == "h1", sex == "t") %>% 
  ggplot()+
  geom_point(aes(age, css, col = method))+
  facet_wrap(~year, scales = "free_y")  

cs3 <- 
  cs2 %>% 
  left_join(pop3) %>% 
  bind_rows(cs %>% mutate(method = "orig")) %>% 
  mutate(inc = 1e5*css/pop) %>% 
  filter(year >= 2009)

cs3 %>% 
  filter(sub == "h3", sex == "t") %>% 
  ggplot()+
  geom_point(aes(age, css, col = method), 
             alpha = .6)+
  facet_wrap(~year, scales = "free_y")  

cs3 %>% 
  filter(sub == "h3", sex == "t", method != "orig") %>% 
  ggplot()+
  geom_point(aes(age, css, col = method), 
             alpha = .6)+
  facet_wrap(~year, scales = "free_y")  +
  labs(title = "h3")

cs3 %>% 
  filter(sub == "h1", sex == "t", method != "orig") %>% 
  ggplot()+
  geom_point(aes(age, css, col = method), 
             alpha = .6)+
  facet_wrap(~year, scales = "free_y")  +
  labs(title = "h1")

# comparing cases between ungrouped and original age groups
comp <- 
  cs2 %>% 
  mutate(age = age-age%%5,
         age = ifelse(age > 85, 85, age)) %>% 
  group_by(year, sex, sub, age, method) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  mutate(css = round(css)) %>% 
  bind_rows(cs %>% 
              mutate(method = "orig") %>% 
              select(-pop) %>% 
              filter(year >= 2009))

comp %>% 
  filter(sub == "h3", sex == "t") %>% 
  ggplot()+
  geom_jitter(aes(age, css, col = method), 
              alpha = .6, height = 0)+
  facet_wrap(~year, scales = "free_y")  +
  theme_bw()+
  labs(title = "h3")

comp %>% 
  filter(sub == "h1", sex == "t") %>% 
  ggplot()+
  geom_jitter(aes(age, css, col = method), 
              alpha = .6, height = 0)+
  facet_wrap(~year, scales = "free_y")  +
  theme_bw()+
  labs(title = "h1")

# weird ungrouping of h1 in 2011
cs3 %>% 
  filter(sub == "h1", sex == "t", year == 2011) %>% 
  ggplot()+
  geom_point(aes(age, css, col = method))


write_rds(cs3, "data_inter/ungrouped_subtypes_2009_2021.rds")
