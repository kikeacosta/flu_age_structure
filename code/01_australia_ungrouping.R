rm(list=ls())
library(tidyverse)
library(ungroup)
options(scipen=999)

dt <- read_csv("data_input/Australie_ka.csv")
pop <- read_csv("data_input/Australia_exposure_hmd.csv")


# population for ungrouping
# only total sex, since 2000 and closing at age 100+
pop2 <- 
  pop %>% 
  filter(Year >= 2000) %>% 
  select(year = Year, age = Age, pop = Total) %>% 
  mutate(age = ifelse(age == "110+", 110, age %>% as.integer()),
         age  = ifelse(age > 100, 100, age)) %>% 
  group_by(year, age) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()

# ages as integers
dt2 <- 
  dt %>% 
  gather(-year, -age, key = subtype, value = cases) %>% 
  separate(age, c("age", "trash")) %>% 
  select(-trash) %>% 
  mutate(age = age %>% as.double(),
         cases = ifelse(is.na(cases) | cases < 0, 0.1, cases))


# re-scaling for flu sub-type imputation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rescale_subtypes <- 
  function(chunk){
    tot <- chunk %>% filter(subtype == "tot") %>% pull(cases)
    chunk %>% 
      filter(subtype != "tot") %>% 
      mutate(cases = (cases/sum(cases))*tot)
  }

# estimating totals for imputation
dt3 <- 
  dt2 %>% 
  group_by(year, age) %>% 
  summarise(cases = sum(cases)) %>% 
  ungroup() %>% 
  mutate(subtype = "tot") %>% 
  bind_rows(dt2) %>% 
  filter(subtype != "Nb A")

# flu sub-type imputation
dt4 <- 
  dt3 %>% 
  group_by(year, age) %>% 
  do(rescale_subtypes(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(subtype = recode(subtype,
                          "Nb H1N1" = "h1",
                          "Nb H3N2" = "h3"))

unique(dt4$subtype)

write_rds(dt4, "data_inter/subtypes_age_australia_std.rds")

# ungrouping ages in single-year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ungroup_this <- 
  function(chunk,
           lambda = 100){
    
    ages <- chunk$age
    # increasing numbers to avoid inconsistent outcomes (decreasing afterwards)
    cass <- chunk$cases*1e5
    # closing at age 100+
    nlast <- 101 - max(ages)
    yr <- unique(chunk$year)
    pops <- 
      pop2 %>% 
      filter(year == yr) %>% 
      pull(pop)
    
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
             cases1 = V1/1e5,
             cases2 = V2/1e5,
             cases3 = V3/1e5,
             cases4 = V4/1e5)
    return(out)
  }

dt5 <- 
  dt4 %>% 
  arrange(year, subtype, age) %>% 
  group_by(year, subtype) %>% 
  do(ungroup_this(chunk = .data, lambda = 1e5)) %>% 
  ungroup()

dt6 <- 
  dt5 %>% 
  gather(starts_with("cases"), key = method, value = cases)

dt6 %>% 
  ggplot()+
  geom_point(aes(age, cases, col = method))+
  facet_grid(subtype~factor(year))+
  theme_bw()
  
dt7 <- 
  dt6 %>% 
  filter(method == "cases3") %>% 
  select(-method)

av_ages <- 
  dt7 %>% 
  mutate(age2 = age-age%%5) %>% 
  group_by(year, subtype, age2) %>%
  summarise(av_age = sum(age*cases)/sum(cases))

dt7 %>%
  filter(year == 2017) %>%
  ggplot(aes(age, cases, col = subtype))+
  geom_point()+
  theme_bw()

dt8 <- 
  dt7 %>% 
  left_join(pop2) %>% 
  mutate(rate = 10000 * cases / pop) %>% 
  group_by(year, subtype) %>% 
  mutate(rate_std = rate / sum(rate)) %>% 
  ungroup()

dt8_gr <- 
  dt8 %>% 
  mutate(age = age - age%%5,
         age = ifelse(age >= 85, 85, age)) %>% 
  group_by(year, subtype, age) %>% 
  summarise(cases = sum(cases),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(rate = 10000 * cases / pop) %>% 
  group_by(year, subtype) %>% 
  mutate(rate_std = rate / sum(rate)) %>% 
  ungroup()

h3_ratio_gr <- 
  dt8_gr %>% 
  select(year, subtype, age, rate_std) %>% 
  spread(subtype, rate_std) %>% 
  mutate(h3_rat = h3 / h1)


h3_ratio_gr %>% 
  ggplot()+
  geom_line(aes(age, h3_rat))+
  scale_y_log10()+
  facet_grid(~year)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

h3_ratio_gr %>% 
  ggplot()+
  geom_line(aes(age, h3_rat, col = factor(year)))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(1910, 2010, 5))+
  # facet_grid(~year)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

h3_ratio <- 
  dt8 %>% 
  select(year, subtype, age, rate_std) %>% 
  spread(subtype, rate_std) %>% 
  mutate(h3_rat = h3 / h1,
         cohort = year - age)

h3_ratio %>% 
  ggplot()+
  geom_line(aes(age, h3_rat))+
  scale_y_log10()+
  facet_grid(~year)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

h3_ratio %>% 
  ggplot()+
  geom_line(aes(cohort, h3_rat, col = factor(year)))+
  scale_y_log10()+
  # facet_grid(~year)+
  scale_x_continuous(breaks = seq(1910, 2010, 5))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

ggsave("figures/h3_ratio_australia_single-year_age.png",
       w = 10, h = 4)
