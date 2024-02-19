rm(list=ls())
library(tidyverse)
library(mgcv)
options(scipen=999)

dt <- read_rds("data_inter/brazil_monthly_baselines.rds")

dts_age <- 
  dt %>% 
  mutate(
    period = case_when(
      # year <= 2008 ~ "sesh3",
      year == 2009 ~ "pan09",
      year == 2012 ~ "wav12",
      year == 2013 ~ "wav13",
      year == 2014 ~ "wav14",
      year == 2015 ~ "wav15",
      year == 2016 ~ "wav16",
      year == 2017 ~ "wav17",
      year == 2018 ~ "wav18",
      year == 2019 ~ "wav19",
      TRUE ~ "seas"
    )
  ) %>% 
  mutate(
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = ifelse(dts > bsn_uc, dts - bsn, 0),
    flu_r = 1e5*flu/exposure
  ) %>% 
  summarise(
    year = mean(year),
    dts = mean(dts),
    bsn = mean(bsn),
    flu = mean(flu), 
    flu_r = mean(flu_r),
    exposure = mean(exposure),
    .by = c(period, cause, sex, age)) %>% 
  mutate(psc = (bsn + flu)/bsn)

# smoothing rates over ages
smooth_age <- function(chunk){
  model <- 
    gam(flu ~ 
          s(age, bs = 'ps', m = c(2,2), k = 20) + 
          offset(log(exposure)), 
        # weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"), gamma = 1e-5)
  
  test <- 
    try(
      res <- 
        predict(model, 
                newdata = chunk,
                type = "response", 
                se.fit = TRUE)
    )
  
  try(
    chunk2 <- 
      chunk %>% 
      mutate(flu_smt = res$fit,
             flu_smt_lc = bsn - 1.96 * res$se.fit,
             flu_smt_uc = bsn + 1.96 * res$se.fit)
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(flu_smt = NA,
             flu_smt_lc = NA,
             flu_smt_uc = NA)
  }
  return(chunk2)
}

wvs <- c("seas", "pan09", "wav13", "wav16", "wav17")

dts_age2 <- 
  dts_age %>% 
  filter(period %in% wvs) %>% 
  mutate(period = factor(period, 
                         levels = c("seas", "pan09", "wav13", "wav16", "wav17"))) %>% 
  group_by(period, cause, sex) %>% 
  do(smooth_age(chunk = .data)) %>% 
  mutate(flu_smt_r = 1e5*flu_smt/exposure) %>% 
  ungroup()

cols <- 
  c("seas" = "#386641",
    "pan09" = "#6a040f",
    "wav13" = "#d00000", 
    "wav16" = "#e85d04", 
    "wav17" = "#3a86ff")

# relative risks by wave subtype ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# seasonal deaths
seas <- 
  dts_age2 %>% 
  filter(period == "seas") %>% 
  select(cause, sex, age, ref_flu_smt_r = flu_smt_r, ref_flu_r = flu_r)

dts_age3 <- 
  dts_age2 %>% 
  filter(period != "seas") %>% 
  select(period, year, cause, sex, age, flu_r, flu_smt_r) %>% 
  left_join(seas) %>% 
  mutate(rr = flu_smt_r / ref_flu_smt_r,
         rr_raw = flu_r / ref_flu_r,
         cohort = year - age) %>% 
  filter(
    sex == "t",
    cause == "pi",
    period %in% c("pan09", "wav17"),
    cohort %in% 1925:2010
  ) %>% 
  group_by(period) %>% 
  mutate(rr_st = rr/mean(rr),
         rr_st_log = log(rr_st),
         rr_st_log2 = rr_st_log + abs(min(rr_st_log)),
         rr_st2 = rr - min(rr),
         ) %>% 
  ungroup()


cols <- 
  c("pan09" = "#d00000",
    "wav17" = "#3a86ff")

bks <- c(0.6, 0.8, 1, 1.2, 1.5, 2, 2.5)
lbs <- 100*(bks-1)

dts_age3 %>% 
  # mutate(rr = ) %>% 
  ggplot()+
  geom_line(aes(cohort, rr_st, col = period), size = 1)+
  scale_x_continuous(breaks = seq(1900, 2010, 10))+
  scale_y_log10(breaks = bks, labels = lbs)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols, labels = c("H1N1", "H3N2"))+
  labs(col = "Year", y = "relative risks (%)",
       title = "H1N1 and H3N2 death risks, compared to a 'typical' flu")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12))
  
ggsave("figures/brazil/ced_slide/serf_rr_coh_vacc.png",
       w = 8, h = 4)

shr <- 
  dts_age3 %>% 
  select(period, cohort, rr_st_log2) %>% 
  spread(period, rr_st_log2) %>% 
  replace_na(list(pan09 = 0,
                  wav17 = 0)) %>% 
  mutate(h1 = pan09/(pan09+wav17),
         h3 = wav17/(pan09+wav17)) %>% 
  select(-pan09, -wav17) %>% 
  gather(h1, h3, key = sub, value = share)
  

shr %>% 
  ggplot()+ 
  geom_area(aes(cohort, share, fill = sub))+
  scale_x_continuous(breaks = seq(1900, 2010, 10))+
  scale_fill_manual(values = c("#d00000", "#3a86ff"))+
  coord_cartesian(expand = 0)+
  labs(fill = "Flu subtype", y = "H1/H3 share")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

ggsave("figures/brazil/ced_slide/serf_rr_coh_vacc_shr.png",
       w = 8, h = 4)

shr %>% 
  ggplot()+ 
  geom_area(aes(cohort, share, fill = sub))+
  scale_x_continuous(breaks = seq(1900, 2010, 10))+
  scale_fill_manual(values = c("#d00000", "#3a86ff"))+
  coord_cartesian(expand = 0)+
  labs(fill = "Flu subtype", y = "H1/H3 share")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12))

ggsave("figures/brazil/ced_slide/legend.png",
       w = 4, h = 2)

