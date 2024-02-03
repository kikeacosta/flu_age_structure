rm(list=ls())
library(tidyverse)
library(mgcv)
options(scipen=999)

dt <- read_rds("data_inter/brazil_monthly_master.rds")

dt2 <- 
  dt %>% 
  filter(year < 2020) %>% 
  mutate(t = 1:n(), 
         w = ifelse(month %in% 4:9, 0, 1),
         sn2 = sin((2*pi*t)/12),
         cs2 = cos((2*pi*t)/12),
         sn4 = sin((4*pi*t)/12),
         cs4 = cos((4*pi*t)/12),
         sn8 = sin((8*pi*t)/12),
         cs8 = cos((8*pi*t)/12),
         sn10 = sin((10*pi*t)/12),
         cs10 = cos((10*pi*t)/12),
         .by = c(cause, sex, age)) %>% 
  rename(exposure = pop/12)

chunk <- 
  dt2 %>% 
  filter(age == 70, sex == "t", cause == "pi")

mod_epi <- 
  gam(dts ~ 
        s(t, bs = 'ps', m = c(2,2)) + 
        sn2 + cs2 +
        sn4 + cs4 +
        sn8 + cs8 +
        sn10 + cs10 +
        offset(log(exposure)), 
      weights = w,
      data = chunk, 
      family = quasipoisson(link = "log"))

res1 <- 
  predict(mod_epi, 
          newdata = chunk,
          type = "response", 
          se.fit = TRUE)

out_1 <- 
  chunk %>% 
  mutate(
    bsn1 = res1$fit,
    bsn1_lc = bsn1 - 1.96 * res1$se.fit,
    bsn1_uc = bsn1 + 1.96 * res1$se.fit
    )

chunk2 <- 
  out_1 %>% 
  mutate(w = ifelse(dts > bsn1_uc, 0, 1))

mod_epi2 <- 
  gam(dts ~ 
        s(t, bs = 'ps', m = c(2,2)) + 
        sn2 + cs2 +
        sn4 + cs4 +
        sn8 + cs8 +
        sn10 + cs10 +
        # s(month, bs = 'cp') +
        offset(log(exposure)), 
      weights = w,
      data = chunk2, 
      family = quasipoisson(link = "log"))

res2 <- 
  predict(mod_epi2, 
          newdata = chunk2,
          type = "response", 
          se.fit = TRUE)

out_2 <- 
  chunk2 %>% 
  mutate(
    bsn2 = res2$fit,
    bsn2_lc = bsn2 - 1.96 * res2$se.fit,
    bsn2_uc = bsn2 + 1.96 * res2$se.fit
  )

out_2 %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = bsn1_lc, ymax = bsn1_uc), 
              alpha = 0.2, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn1), col = "red")+
  geom_line(aes(date, bsn2), col = "blue")+
  theme_bw()



dt3 <- 
  dt2 %>% 
  group_by(cause, sex, age) %>% 
  do(est_epi_periods(chunk = .data)) %>% 
  ungroup()

dt4 <- 
  dt3 %>% 
  mutate(epi = ifelse(dts > bsn_uc, 1, 0))

dt4 %>% 
  filter(sex == "t",
         age == 80) %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = bsn_lc, ymax = bsn_uc), 
              alpha = 0.2, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "red")+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

temp <- dt4

temp2 <- 
  temp %>% 
  select(year, month, date, cause, sex, age, bsn_1 = bsn) %>% 
  left_join(dt4 %>% 
              select(year, month, date, cause, sex, age, bsn_2 = bsn))
              


dt3 <- 
  dt2 %>% 
  group_by(cause, sex, age) %>% 
  do(est_mth_baseline_pi(chunk = .data)) %>% 
  ungroup()

# write_rds(dt3, "data_inter/brazil_monthly_baselines.rds",
#           compress = "xz")


# functions for baseline estimation ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
est_bsns <- function(chunk){
  
  mod_sin <- 
    gam(dts ~ 
          s(t, bs = 'ps', m = c(2,2)) + 
          sn2 + cs2 +
          offset(log(exposure)), 
        weights = w, 
        data = chunk, 
        family = quasipoisson(link = "log"))
  
  mod_spl <- 
    gam(dts ~ 
          s(t, bs = 'ps', m = c(2,2)) + 
          s(month, bs = 'cp') +
          offset(log(exposure)), 
        weights = w, 
        data = chunk, 
        family = quasipoisson(link = "log"))
  
  test_sin <- 
    try(
      res_sin <- 
        predict(mod_sin, 
                newdata = chunk %>% mutate(w = 1),
                type = "response", 
                se.fit = TRUE)
    )
  test_spl <- 
    try(
      res_spl <- 
        predict(mod_spl, 
                newdata = chunk %>% mutate(w = 1),
                type = "response", 
                se.fit = TRUE)
    )
  
  try(
    chunk2 <- 
      chunk %>% 
      mutate(bsn_sin = res_sin$fit,
             bsn_spl = res_spl$fit)
  )
  
  if(class(test_sin) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(bsn_sin = NA)
  }
  if(class(test_spl) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(bsn_spl = NA)
  }
  return(chunk2)
}
  



est_mth_baseline_pi <- function(chunk){
  
  model <- 
    gam(dts ~ 
          s(t, bs = 'ps', m = c(2,2)) + 
          sn2 + cs2 +
          sn4 + cs4 +
          sn8 + cs8 +
          sn10 + cs10 +
          offset(log(exposure)), 
        weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"))
  
  test <- 
    try(
      res <- 
        predict(model, 
                newdata = chunk %>% mutate(w = 1),
                type = "response", 
                se.fit = TRUE)
    )
  
  try(
    chunk2 <- 
      chunk %>% 
      mutate(bsn = res$fit,
             bsn_lc = bsn - 1.96 * res$se.fit,
             bsn_uc = bsn + 1.96 * res$se.fit)
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(bsn = NA,
             bsn_lc = NA,
             bsn_uc = NA)
  }
  
  return(chunk2)
}
