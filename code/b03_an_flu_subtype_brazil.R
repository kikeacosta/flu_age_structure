rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
options(scipen=999)

flu <- 
  read_rds("data_inter/flu_data_brazil_exposures_2009_2019.rds")



flu %>% 
  group_by(sub) %>% 
  summarise(sum(css))

flu %>% 
  group_by(year, sub) %>% 
  summarise(css = sum(css)) %>% 
  spread(sub, css)




# H1 ====
# ~~~~~~~
h1 <- 
  flu %>% 
  filter(sub == "h1",
         sex == "t") %>% 
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
  flu %>% 
  filter(sub == "h3",
         sex == "t") %>% 
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
  flu %>% 
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

