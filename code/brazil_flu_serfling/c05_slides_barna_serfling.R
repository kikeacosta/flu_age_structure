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
        family = quasipoisson(link = "log"), gamma = 1e-2)
  
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


# death rates by season and waves ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts_age2 %>%
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "seas", "wav13", "wav16", "wav17")) %>%
  ggplot()+
  geom_point(aes(age, flu_r, col = period), alpha = 0.4, size = 0.5)+
  geom_line(aes(age, flu_smt_r, col = period))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "death rates (/100K)",
       title = "H1 and H3 seasons")+
  # facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/bcn_slides/serf_mx_h1_h3_age.png",
       w = 8, h = 4)

dts_age2 %>%
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "seas", "wav13", "wav16")) %>%
  ggplot()+
  geom_point(aes(age, flu_r, col = period), alpha = 0.4, size = 0.5)+
  geom_line(aes(age, flu_smt_r, col = period))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "death rates (/100K)",
       title = "H1p waves")+
  # facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/bcn_slides/serf_mx_h1_age.png",
       w = 8, h = 4)


dts_age2 %>%
  filter(sex == "t",
         cause == "pi",
         period %in% c("seas", "wav17")) %>%
  ggplot()+
  geom_point(aes(age, flu_r, col = period), alpha = 0.4, size = 0.5)+
  geom_line(aes(age, flu_smt_r, col = period))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "death rates (/100K)",
       title = "H3 season")+
  # facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/bcn_slides/serf_mx_h3_age.png",
       w = 8, h = 4)



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
         cohort = year - age)

dts_age3 %>% 
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "seas", "wav13", "wav16"),
         cohort %in% 1925:2010) %>% 
  ggplot()+
  geom_line(aes(age, rr, col = period))+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "relative risks",
       title = "H1 waves")+
  theme_bw()

ggsave("figures/brazil/bcn_slides/serf_rr_h1_age.png",
       w = 8, h = 4)

h1_coh <- 
  dts_age3 %>% 
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "seas", "wav13", "wav16"),
         cohort %in% 1925:2010) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "relative risks",
       title = "H1 waves")+
  theme_bw();h1_coh

ggsave("figures/brazil/bcn_slides/serf_rr_h1_coh.png",
       w = 8, h = 4)


h1_coh+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/serf_rr_h1_coh_lines.png",
       w = 8, h = 4)



dts_age3 %>% 
  filter(sex == "t",
         cause == "pi",
         period %in% c("wav17"),
         cohort %in% 1925:2010) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "relative risks",
       title = "H3 wave relative risks")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/serf_rr_h3_coh.png",
       w = 8, h = 4)


dts_age3 %>% 
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "wav17"),
         cohort %in% 1925:2010) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "relative risks",
       title = "H1 pandemic vs. H3 wave")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/serf_rr_h1_h3_coh.png",
       w = 8, h = 4)



# testing mirror
dts_age3 %>% 
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "wav17"),
         cohort %in% 1925:2010) %>% 
  mutate(rr = ifelse(period == "wav17", 1/rr, rr)) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "relative risks",
       title = "H1 pandemic vs. H3 wave")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/serf_rr_h1_h3_coh_test.png",
       w = 8, h = 4)

