rm(list=ls())
source("code/brazil_flu_incidence/b00_functions.R")

# exposures
data(popAge1dt)

pop <- 
  popAge1dt %>% 
  select(year, name, age, m = popM, f = popF, t = pop) %>% 
  filter(name == "Brazil",
         year %in% 1999:2022) %>% 
  pivot_longer(c(m, f, t), 
               names_to = "sex", values_to = "pop") %>% 
  select(year, age, sex, pop) %>% 
  mutate(pop = pop*1e3)

unique(pop$sex)

# ungrouping exposures at ages 90-100
chunk <- 
  pop %>% 
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

pop2 <- 
  pop %>% 
  group_by(year, sex) %>% 
  do(ungroup_open_age(chunk = .data)) %>% 
  ungroup() 

write_rds(pop2, 
          "data_inter/brazil_exposures_1999_2022_ages_0_110.rds",
          compress = "xz")
