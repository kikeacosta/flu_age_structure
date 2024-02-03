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
      year == 2016 ~ "wav16",
      year == 2018 ~ "wav18",
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

dts_age %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(age, flu_r, col = period))+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_age %>%
  filter(sex == "t") %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point()+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  geom_smooth(span = 0.1, alpha = 0.1)+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_age %>%
  filter(sex == "t",
         cause == "pi") %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point()+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  geom_smooth(span = 0.1, alpha = 0.1)+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_age %>%
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "oth", "wav13")) %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point()+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  geom_smooth(span = 0.1, alpha = 0.1)+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# smoothing rates over ages
smooth_age <- function(chunk){
  model <- 
    gam(flu ~ 
          s(age, bs = 'ps', m = c(2,2), k = 30) + 
          offset(log(exposure)), 
        # weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"), gamma = 1e-10)
  
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

dts_age2 <- 
  dts_age %>% 
  group_by(period, cause, sex) %>% 
  do(smooth_age(chunk = .data)) %>% 
  mutate(flu_smt_r = 1e5*flu_smt/exposure) %>% 
  ungroup()

dts_age2 %>%
  filter(sex == "t",
         cause == "pi",
         period %in% c("pan09", "seas", "wav13")) %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point(size = .3)+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  # geom_smooth(span = 0.1, alpha = 0.1)+
  geom_line(aes(age, flu_smt_r, col = period))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/age_structure_pi.png")

dts_age2 %>%
  filter(sex == "t",
         cause %in% c("pi", "cvd", "res"),
         period %in% c("pan09", "seas", "wav13")) %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point(size = .15, alpha = 0.8)+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  # geom_smooth(span = 0.1, alpha = 0.1)+
  geom_line(aes(age, flu_smt_r, col = period), alpha = 0.4)+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/age_structure.png",
       w = 10, 
       h = 4)

dts_age2 %>%
  filter(
    sex == "t",
    cause %in% c("pi", "cvd", "res"),
    # period %in% c("pan09", "seas", "wav13")
    ) %>%
  ggplot(aes(age, flu_r, col = period))+
  geom_point(size = .15, alpha = 0.8)+
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 8), se = FALSE)+
  # geom_smooth(span = 0.1, alpha = 0.1)+
  geom_line(aes(age, flu_smt_r, col = period), alpha = 0.4)+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

ggsave("figures/brazil/age_structure_waves.png",
       w = 10, 
       h = 4)

seas <- 
  dts_age2 %>% 
  filter(period == "seas") %>% 
  select(cause, sex, age, ref_flu_smt_r = flu_smt_r)

dts_age3 <- 
  dts_age2 %>% 
  filter(period != "seas") %>% 
  select(period, year, cause, sex, age, flu_smt_r) %>% 
  left_join(seas) %>% 
  mutate(rr = flu_smt_r / ref_flu_smt_r,
         cohort = year - age)

dts_age3 %>% 
  filter(sex == "t") %>% 
  ggplot()+
  geom_line(aes(age, rr, col = period))+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_age3 %>% 
  filter(sex == "t",
         period %in% c("pan09", "wav13")) %>% 
  ggplot()+
  geom_line(aes(age, rr, col = period))+
  scale_x_continuous(breaks = seq(0,100, 10))+
  scale_y_log10()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~cause)+
  theme_bw()

dts_age3 %>% 
  filter(sex == "t",
         period %in% c("pan09", "wav13", "wav16")) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900,2010, 10))+
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = c(1918, 1957, 1968, 1977), linetype = "dashed")+
  facet_wrap(~cause)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7))

ggsave("figures/brazil/cohort_rr.png",
       w = 10, 
       h = 5)

dts_age3 %>% 
  filter(
    sex == "t",
    # period %in% c("pan09", "wav13", "wav16")
  ) %>% 
  ggplot()+
  geom_line(aes(cohort, rr, col = period))+
  scale_x_reverse(breaks = seq(1900,2010, 10))+
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = c(1918, 1957, 1968, 1977), linetype = "dashed")+
  facet_wrap(~cause)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7),
        legend.position = "bottom")

ggsave("figures/brazil/cohort_rr_waves.png",
       w = 10, 
       h = 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




ag <- 90

# dt %>%
#   filter(age == ag, sex == "t") %>% 
#   ggplot()+
#   geom_ribbon(aes(date, ymin = bsn_lp, ymax = bsn_up), 
#               alpha = 0.2, fill = "red")+
#   geom_line(aes(date, dts))+
#   geom_line(aes(date, bsn), col = "red")+
#   labs(title = ag)+
#   facet_wrap(~cause, scales = "free_y")+
#   theme_bw()

# identifying months of flu activity during the pandemic
dt %>%
  filter(age == ag, sex == "t",
         year == 2013) %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = bsn_lp, ymax = bsn_up), 
              alpha = 0.2, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "red")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(title = ag)+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# june - september 2009
int09 <- interval(ymd("2009-06-01"), ymd("2009-10-31"))

# identifying months of flu activity during the pandemic
dt %>%
  filter(age == ag, sex == "t",
         year == 2013) %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = bsn_lc, ymax = bsn_uc), 
              alpha = 0.2, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "red")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  labs(title = ag)+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

int09 <- interval(ymd("2009-06-01"), ymd("2009-10-31"))
# june - september 2009


# identifying months of flu activity during before H1 pandemic
dt %>%
  filter(age == ag, sex == "t",
         year %in% 2000:2008) %>% 
  ggplot()+
  # geom_ribbon(aes(date, ymin = bsn_lp, ymax = bsn_up), 
  #             alpha = 0.2, fill = "red")+
  geom_line(aes(month, dts, col = factor(year), group = year))+
  # geom_line(aes(date, bsn), col = "red")+
  scale_x_continuous(breaks = 1:12)+
  labs(title = ag)+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# identifying months of flu activity during H1 pandemic and waves
dt %>%
  filter(age == ag, sex == "t",
         year %in% c(2009, 2013, 2016, 2018)) %>% 
  ggplot()+
  # geom_ribbon(aes(date, ymin = bsn_lp, ymax = bsn_up), 
  #             alpha = 0.2, fill = "red")+
  geom_line(aes(month, dts, col = factor(year), group = year))+
  # geom_line(aes(date, bsn), col = "red")+
  scale_x_continuous(breaks = 1:12)+
  labs(title = ag)+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# monthly deaths during 2009 pandemic
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inth3 <- interval(ymd("2000-04-01"), ymd("2008-09-30"))
int13 <- interval(ymd("2013-04-01"), ymd("2013-08-31"))
int16 <- interval(ymd("2016-03-01"), ymd("2016-06-30"))
int18 <- interval(ymd("2018-05-01"), ymd("2018-08-31"))

dts_age <- 
  dt %>% 
  mutate(period = case_when(date %within% inth3 ~ "sesh3",
                            date %within% int09 ~ "pan09",
                            date %within% int13 ~ "wav13",
                            date %within% int16 ~ "wav16",
                            date %within% int18 ~ "wav18",
                            TRUE ~ "oth")) %>% 
  mutate(
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = ifelse(dts > bsn_up, dts - bsn, 0),
    flu_r = 1e5*flu/exposure
  ) %>% 
  summarise(dts = mean(dts),
            bsn = mean(bsn),
            flu = mean(flu), 
            flu_r = mean(flu_r),
            exposure = mean(exposure),
            .by = c(period, cause, sex, age))

dts_coh <- 
  dt %>% 
  mutate(period = case_when(date %within% inth3 ~ "sesh3",
                            date %within% int09 ~ "pan09",
                            date %within% int13 ~ "wav13",
                            date %within% int16 ~ "wav16",
                            date %within% int18 ~ "wav18",
                            TRUE ~ "oth"),
         cohort = year - age) %>% 
  mutate(
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = ifelse(dts > bsn_up, dts - bsn, 0),
    flu_r = 1e5*flu/exposure
  ) %>% 
  summarise(dts = mean(dts),
            bsn = mean(bsn),
            flu = mean(flu), 
            flu_r = mean(flu_r),
            exposure = mean(exposure),
            .by = c(period, cause, sex, cohort))

dts_age %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(age, flu_r, col = period))+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_coh %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(cohort, flu_r, col = period))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_coh %>%
  filter(sex == "t",
         period == "sesh3") %>%
  ggplot()+
  geom_line(aes(cohort, flu_r, col = period))+
  geom_vline(xintercept = c(1957, 1968, 1977), linetype = "dashed")+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(1910, 2010, 10))+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

dts_age %>%
  filter(sex == "t",
         period %in% c("sesh3", "pan09")) %>%
  ggplot()+
  geom_line(aes(age, flu_r, col = period))+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# dts_p <- 
#   dt %>% 
#   filter(year == 2009,
#          month %in% 7:10) %>% 
#   mutate(
#     bsn = ifelse(bsn > dts, dts, bsn),
#     flu = ifelse(dts > bsn_up, dts - bsn, 0),
#     flu_r = 1e5*flu/exposure
#   ) %>% 
#   summarise(dts = mean(dts),
#             bsn = mean(bsn),
#             flu = mean(flu), 
#             flu_r = mean(flu_r),
#             exposure = mean(exposure),
#             .by = c(cause, sex, age)) %>% 
#   mutate(type = "pandemic 2009")
# 
# dts_13w <- 
#   dt %>% 
#   filter(year == 2013,
#          month %in% 4:8) %>% 
#   mutate(
#     bsn = ifelse(bsn > dts, dts, bsn),
#     flu = ifelse(dts > bsn_up, dts - bsn, 0),
#     flu_r = 1e5*flu/exposure
#   ) %>% 
#   summarise(dts = mean(dts),
#             bsn = mean(bsn),
#             flu = mean(flu), 
#             flu_r = mean(flu_r),
#             exposure = mean(exposure),
#             .by = c(cause, sex, age)) %>% 
#   mutate(type = "wave 2013")
# 
# dts_16w <- 
#   dt %>% 
#   filter(year == 2016,
#          month %in% 2:8) %>% 
#   mutate(
#     bsn = ifelse(bsn > dts, dts, bsn),
#     flu = ifelse(dts > bsn_up, dts - bsn, 0),
#     flu_r = 1e5*flu/exposure
#   ) %>% 
#   summarise(dts = mean(dts),
#             bsn = mean(bsn),
#             flu = mean(flu), 
#             flu_r = mean(flu_r),
#             exposure = mean(exposure),
#             .by = c(cause, sex, age)) %>% 
#   mutate(type = "wave 2016")

# monthly deaths during seasonal outbreaks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts_s <- 
  dt %>% 
  filter(
    # year != 2009,
    month %in% 4:9,
    year %in% 2000:2008
  ) %>% 
  mutate(
    bsn = ifelse(bsn > dts, dts, bsn),
    flu = ifelse(dts > bsn_up, dts - bsn, 0),
    flu_r = 1e5*flu/exposure
  )%>% 
  summarise(dts = mean(dts),
            bsn = mean(bsn),
            flu = mean(flu), 
            flu_r = mean(flu_r),
            exposure = mean(exposure),
            .by = c(cause, sex, age)) %>% 
  mutate(type = "seasonal average")


dts <- 
  bind_rows(dts_p, dts_s, dts_13w, dts_16w)

# plotting pandemic vs seasonal mortality
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dts %>%
  filter(sex == "t") %>%
  ggplot()+
  geom_line(aes(age, flu_r, col = type))+
  scale_y_log10()+
  facet_wrap(~cause, scales = "free_y")+
  theme_bw()

# plotting p-scores
# ~~~~~~~~~~~~~~~~~~~
flu_p %>%
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


