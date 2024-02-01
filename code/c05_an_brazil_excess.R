rm(list=ls())
library(tidyverse)
# library(mgcv)
options(scipen=999)

dt <- read_rds("data_inter/brazil_monthly_baselines.rds")


ag <- 52
dt %>%
  filter(age == ag, sex == "t") %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = bsn_lp, ymax = bsn_up), alpha = 0.2, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "red")+
  # geom_line(aes(date, bsn_sin), col = "blue")+
  labs(title = ag)+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

pn1 <- 
  dt %>% 
  filter(year == 2009,
         month %in% 9:12) %>% 
  mutate(
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = dts - bsn) %>% 
  summarise(dts = sum(dts),
            bsn = sum(bsn),
            flu = sum(flu), 
            exposure = sum(exposure),
            .by = c(cause, sex, age)) %>% 
  mutate(psc = dts/bsn)

pn1 %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(age, flu))+
  # geom_line(aes(date, bsn), col = "red")+
  # geom_line(aes(date, bsn_sin), col = "blue")+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

pn1 %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(age, psc))+
  # geom_line(aes(date, bsn), col = "red")+
  # geom_line(aes(date, bsn_sin), col = "blue")+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()




pn1 <- 
  dt %>% 
  # filter(year == 2009,
  #        month %in% 9:12) %>% 
  mutate(
    evt = case_when()
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = dts - bsn) %>% 
  summarise(dts = sum(dts),
            bsn = sum(bsn),
            flu = sum(flu), 
            exposure = sum(exposure),
            .by = c(cause, sex, age)) %>% 
  mutate(psc = dts/bsn)


