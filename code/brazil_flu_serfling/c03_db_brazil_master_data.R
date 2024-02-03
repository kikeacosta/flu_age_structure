rm(list=ls())
library(tidyverse)
options(scipen=999)

dts <- read_rds("data_inter/brazil_monthly_cause_death.rds")
pop <- read_rds("data_inter/brazil_monthly_pop.rds")

# adding cardio_resp
dts2 <- 
  dts %>% 
  filter(cause %in% c("cvd", "res", "pi")) %>% 
  reframe(dts = sum(dts),
          .by = date, year, month, sex, age) %>% 
  mutate(cause = "cvd_res") %>% 
  bind_rows(dts)

dt <- 
  dts2 %>% 
  left_join(pop) %>% 
  mutate(mx = dts*1e5/pop)


# visualizing annual values and weekly estimates
dt %>% 
  filter(age == 52) %>% 
  ggplot()+
  geom_line(aes(date, mx), linewidth = 1)+
  # geom_point(aes(date, pop), col = "red", size = 1.5, alpha = 0.7)+
  facet_grid(cause~sex, scales = "free_y")+
  theme_bw()

write_rds(dt, "data_inter/brazil_monthly_master.rds",
          compress = "xz")
