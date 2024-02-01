rm(list=ls())
library(tidyverse)
library(wpp2022)

# exposures
data(popAge1dt)

pop <- 
  popAge1dt %>% 
  select(year, name, age, m = popM, f = popF, t = pop) %>% 
  filter(name == "Brazil",
         year %in% 1999:2020) %>% 
  pivot_longer(c(m, f, t), 
               names_to = "sex", values_to = "pop") %>% 
  select(year, age, sex, pop) %>% 
  mutate(month = 6)

unique(pop$age)

# interpolating population into months


# from annual to weekly exposures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wee need weekly exposures between 2010 and 2023
# then, we can interpolate between 2009 and 2024

# first, create a grid with all weeks that we need to populate with 
# weekly exposures
pop_empty <- 
  # 52 weeks per year
  expand_grid(year = 1999:2020, month = 1:12, sex = c("m", "f", "t"), age = 0:100) %>%
  mutate(date = make_date(d = 15, m = month, y = year)) %>% 
  # arrange it chronologically
  arrange(date, sex, age) 
  
# preparing data for interpolation
pop_to_int <- 
  pop_empty %>%  
  # adding annual population in midyear to week 26 (the midyear week!)
  left_join(pop) %>% 
  mutate(t = 1:n(), .by = c(sex, age))

inter_pop <- function(dt){

  # extracting values for interpolation
  xs <- dt %>% drop_na() %>% pull(t) # weeks with data
  ys <- dt %>% drop_na() %>% pull(pop) # available pop data
  ts <- dt %>% pull(t) # all weeks in which we need estimates
  
  # smoothing using cubic splines
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # the "spline()" function allows to interpolate by constraining the curve 
  # to match the original data. In other words, it is strictly interpolation 
  # rather than smoothing  
  
  # extracting predictions from the model
  dt2 <- 
    dt %>% 
    mutate(pop_int = spline(xs, ys, xout = ts)$y)

  return(dt2)
}

# weekly exposures to use
pop_int <- 
  pop_to_int %>% 
  group_by(sex, age) %>% 
  do(inter_pop(dt = .data)) %>% 
  ungroup()


# visualizing annual values and weekly estimates
pop_int %>% 
  filter(age == 20) %>% 
  ggplot()+
  geom_line(aes(date, pop_int), linewidth = 1)+
  geom_point(aes(date, pop), col = "red", size = 1.5, alpha = 0.7)+
  facet_grid(~sex)+
  theme_bw()

pop_out <- 
    pop_int %>% 
    select(-pop, -t) %>% 
    rename(pop = pop_int)
  
write_rds(pop_out, "data_inter/brazil_monthly_pop.rds",
          compress = "xz")
