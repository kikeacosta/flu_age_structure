rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

flu <- 
  read_rds("data_inter/flu_data_brazil_exposures_2009_2019.rds")

flu %>% 
  group_by(sub, outcome) %>% 
  summarise(sum(value))

flu %>% 
  group_by(year, sub) %>% 
  summarise(css = sum(css)) %>% 
  spread(sub, css)

# H1 ====
# ~~~~~~~
h1 <- 
  flu %>% 
  filter(sub == "h1",
         sex == "t",
         outcome == "cases") %>% 
  group_by(cohort, sex) %>% 
  summarise(value = sum(value),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(mx = 1e5 * value / pop)

h1 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, value))+
  geom_line(aes(cohort, value))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()+
  labs(title = "H1")

h1 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, mx))+
  geom_line(aes(cohort, mx))+
  geom_vline(xintercept = c(1918, 1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "H1")


# H3 ====
# ~~~~~~~
h3 <- 
  flu %>% 
  filter(sub == "h3",
         sex == "t",
         outcome == "cases") %>% 
  group_by(cohort, sex) %>% 
  summarise(value = sum(value),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(mx = 1e5 * value / pop)

h3 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, value))+
  geom_line(aes(cohort, value))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  theme_bw()+
  labs(title = "H3")

h3 %>% 
  filter(cohort %in% 1910:2012) %>% 
  ggplot()+
  geom_point(aes(cohort, mx))+
  geom_line(aes(cohort, mx))+
  geom_vline(xintercept = c(1918, 1957, 1968, 2009), linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "H3")

hs <- 
  flu %>% 
  filter(sub %in% c("h1", "h3"), 
         sex == "t",
         outcome == "cases") %>% 
  group_by(sub, cohort, outcome) %>% 
  summarise(value = sum(value),
            pop = sum(pop)) %>% 
  ungroup() %>% 
  mutate(mx = 1e5*value/pop) %>% 
  group_by(sub, outcome) %>% 
  mutate(cx_cs = value / sum(value),
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

rh1h3 <- 
  hs %>% 
  select(sub, cohort, cx_cs) %>% 
  spread(sub, cx_cs) %>% 
  replace_na(list(h1 = 0, h3 = 0)) %>% 
  mutate(rh1h3 = h1/h3) %>% 
  filter(cohort %in% 1915:2012) 


rh1h3 %>% 
  ggplot()+
  geom_point(aes(cohort, rh1h3))+
  geom_line(aes(cohort, rh1h3))+
  geom_vline(xintercept = c(1957, 1968, 2009), linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))+
  scale_y_log10()+
  theme_bw()+
  labs(title = "ratio H1/h3")



