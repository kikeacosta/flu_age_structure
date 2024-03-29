rm(list=ls())
source("code/brazil_flu_incidence/b00_functions.R")

pop <- 
  read_rds("data_inter/brazil_exposures_1999_2022_ages_0_110.rds")

dt <- 
  read_rds("data_inter/flu_data_brazil_2009_2019_v3.rds") %>% 
  mutate(age = year - cohort)

unique(dt$sex)

unique(dt$year)
unique(dt$typ)
table(dt$outcome)
unique(dt$sub)
table(dt$sub)
table(dt$year, dt$sub, dt$outcome)
table(dt$typ, dt$sub)
table(dt$year, dt$sub)
table(dt$hosp)
table(dt$flu, dt$hosp)
table(dt$outcome)
table(dt$year)
table(dt$age)
table(dt$year, dt$sub)
# table(dt2$sub, dt2$outcome)

# keeping only subtype and redistributing deaths by sex
dt2 <- 
  dt %>% 
  mutate(age = year - cohort,
         age = case_when(age < 0 ~ 0, 
                         age > 110 ~ 110,
                         TRUE ~ age)) %>% 
  select(-typ, -cohort) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age, flu, sub, hosp, outcome)) %>% 
  # imputing sex distribution
  spread(sex, value) %>% 
  replace_na(list(f = 0, m = 0, i = 0)) %>% 
  mutate(
    t = f+m+i,
    dts_sum = f+m,
    f = round(f*t/dts_sum),
    m = round(m*t/dts_sum)
  ) %>% 
  select(-dts_sum, -i) %>% 
  gather(t, f, m, key = sex, value = value) %>% 
  mutate(dth = ifelse(outcome %in% c("death_flu", "death_inv", "death_oth"), 
                      1, 0))

# grouping data 
# ili
ili <- 
  dt2 %>%
  filter(sex == "t") %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "ili")

# sari
sari <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(hosp == 1) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "sari")

# flu cases
flu <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1) %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "flu")

# H1 cases
h1 <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1,
         sub == "h1") %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "h1")

# H3 cases
h3 <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu == 1,
         sub == "h3") %>% 
  mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "h3")

non_flu <- 
  dt2 %>% 
  filter(sex == "t") %>% 
  filter(flu != 1,
         sub == "z") %>% 
  # mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "nonflu")

# putting all together
all <- 
  bind_rows(ili, sari, flu, h1, h3, non_flu) %>% 
  complete(year, sex, age, fill = list(css = 0, hsp = 0, dts = 0)) %>% 
  arrange(type, year, sex, age) %>% 
  left_join(pop) %>% 
  filter(year %in% c(2009, 2012, 2013, 2016, 2017, 2018, 2019),
         sex == "t")


# smoothing rates along ages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# incidence
ix <- 
  all %>% 
  rename(value = css, exposure = pop) %>% 
  group_by(year, sex, type) %>% 
  do(smooth_age(chunk = .data)) %>% 
  ungroup() 

# death rates
mx <- 
  all %>% 
  rename(value = dts, exposure = pop) %>% 
  group_by(year, sex, type) %>% 
  do(smooth_age(chunk = .data)) %>% 
  ungroup()

# hosp risk
hsp <- 
  all %>% 
  rename(value = hsp, exposure = pop) %>% 
  group_by(year, sex, type) %>% 
  do(smooth_age(chunk = .data)) %>% 
  ungroup()

# estimating CFR and smoothed CFR from smoothed rates
# CFR
cfr <- 
  all %>% 
  rename(value = dts, exposure = css) %>% 
  left_join(ix %>% 
              select(year, sex, age, type, exp_smt = val_smt)) %>% 
  left_join(mx %>% 
              select(year, sex, age, type, val_smt = val_smt)) %>% 
  mutate(value_r = value/exposure,
         val_smt_r = val_smt/exp_smt)


all %>% 
  filter(type == "h1",
         year == 2009) %>% 
  mutate(cfr = dts/css) %>% 
  ggplot()+
  geom_point(aes(age, cfr), col = "red")
  

# plots ====
# ~~~~~~~~~~
# Incidence
ix %>% 
  filter(sex == "t",
         type %in% c("h1", "h3")) %>%
  ggplot()+
  geom_point(aes(age, value, col = factor(year)))+
  geom_line(aes(age, val_smt, col = factor(year)))+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  facet_wrap(~type, scales = "free_y")

# over ages
ix %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  # geom_point(aes(age, value_r, col = factor(year)))+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  theme_bw()+
  scale_y_log10()

# over cohorts
ix %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt, col = factor(year)))+
  geom_point(aes(cohort, val_smt, col = factor(year)))+
  scale_x_reverse()+
  theme_bw()+
  scale_y_log10()+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed")

cols <- c("#6a040f", "#d00000", "#e85d04")

ix %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1920:2010) %>%
  ggplot()+
  geom_point(aes(age, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  theme_bw()+
  scale_y_log10()+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "incidence (/100K)")+
  scale_x_continuous(breaks = seq(0, 100, 10))
  # geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/surv_ix_h1_age.png",
       w = 8, h = 4)


ix %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1920:2010) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  geom_point(aes(cohort, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  theme_bw()+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "incidence (/100K)")+
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1920, 2010, 10))
ggsave("figures/brazil/bcn_slides/surv_ix_h1_coh.png",
       w = 8, h = 4)

ix %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1920:2010) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  geom_point(aes(cohort, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  theme_bw()+
  scale_color_manual(values = cols)+
  labs(col = "Year", y = "incidence (/100K)")+
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1920, 2010, 10))+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/surv_ix_h1_coh_lines.png",
       w = 8, h = 4)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hospitalization risk
hsp %>% 
  filter(sex == "t",
         type %in% c("h1", "h3")) %>%
  ggplot()+
  geom_point(aes(age, value, col = factor(year)))+
  geom_line(aes(age, val_smt, col = factor(year)))+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  facet_wrap(~type, scales = "free_y")

# over ages
hsp %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  # geom_point(aes(age, value_r, col = factor(year)))+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  theme_bw()+
  scale_y_log10()

# over cohorts
hsp %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt, col = factor(year)))+
  geom_point(aes(cohort, val_smt, col = factor(year)))+
  scale_x_reverse()+
  theme_bw()+
  scale_y_log10()+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed")


hsp %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h3",
         cohort >= 1920) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  theme_bw()+
  scale_y_log10()+
  scale_x_reverse()+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mortality ====
# ~~~~~~~~~~~~~~
# death counts
mx %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  geom_point(aes(age, value, col = factor(year)))+
  geom_line(aes(age, val_smt, col = factor(year)))+
  # facet_wrap(~year, scales = "free_y")+
  theme_bw()

# death rates
mx %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1925) %>%
  ggplot()+
  # geom_point(aes(age, value))+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  scale_y_log10()+
  # facet_wrap(~year, scales = "free_y")+
  theme_bw()

# death rates over cohorts
mx %>% 
  mutate(cohort = year - age,
         value_r = value / exposure,
         val_smt_r = val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1925) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  scale_y_log10()+
  # facet_wrap(~year, scales = "free_y")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed")



mx %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1925:2010) %>%
  ggplot()+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  geom_point(aes(age, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(col = "Year", y = "death rates (/100K)")+
  theme_bw()
ggsave("figures/brazil/bcn_slides/surv_mx_h1_age.png",
       w = 8, h = 4)


mx %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1925:2010) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  geom_point(aes(cohort, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(col = "Year", y = "death rates (/100K)")+
  theme_bw()
ggsave("figures/brazil/bcn_slides/surv_mx_h1_coh.png",
       w = 8, h = 4)


mx %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2013, 2016),
         cohort %in% 1925:2010) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  geom_point(aes(cohort, value_r, col = factor(year)), alpha = 0.4, size = 0.5)+
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(col = "Year", y = "death rates (/100K)")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1984), linetype = "dashed")

ggsave("figures/brazil/bcn_slides/surv_mx_h1_coh_lines.png",
       w = 8, h = 4)






# case fatality rates ====
# ~~~~~~~~~~~~~~~~~~~~~~~~
# over age
cfr %>% 
  mutate(cohort = year - age) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1910) %>%
  ggplot()+
  geom_point(aes(age, value_r, col = factor(year)))+
  geom_line(aes(age, val_smt_r, col = factor(year)))+
  # facet_wrap(~year, scales = "free_y")+
  theme_bw()

# over cohorts
cfr %>% 
  mutate(cohort = year - age) %>% 
  filter(sex == "t",
         type == "h1",
         cohort >= 1920) %>%
  ggplot()+
  geom_point(aes(cohort, value_r, col = factor(year)))+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  # facet_wrap(~year, scales = "free_y")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1985, 2009), linetype = "dashed")

# selected waves
cfr %>% 
  mutate(cohort = year - age) %>% 
  filter(sex == "t",
         type == "h1",
         year %in% c(2009, 2012, 2013, 2016),
         cohort >= 1920) %>%
  ggplot()+
  geom_point(aes(cohort, value_r, col = factor(year)))+
  geom_line(aes(cohort, val_smt_r, col = factor(year)))+
  # facet_wrap(~year, scales = "free_y")+
  scale_x_reverse()+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968, 1984, 2009), linetype = "dashed")

# relative risks
# H1 relative to all not identified as flu
rrh1 <- 
  h1 %>% 
  select(year, age, css_h1 = css, dts_h1 = dts) %>% 
  left_join()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# total flu circulation ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# cases ====
# ~~~~~~~~~~

css_ili <- 
  dt2 %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

css_flu <- 
  dt2 %>% 
  filter(flu == 1) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

css_h1 <- 
  dt2 %>% 
  filter(sub == "h1") %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

css_h3 <- 
  dt2 %>% 
  filter(sub == "h3") %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

# deaths ====
# ~~~~~~~~~~
dths_names <- c("death_flu","death_oth", "death_inv")
unique(dt2$outcome)

dts_ili <- 
  dt2 %>% 
  filter(outcome %in% dths_names) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age)) %>% 
  spread(age, value)

dts_flu <- 
  dt2 %>% 
  filter(flu == 1,
         outcome %in% dths_names) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

dts_h1 <- 
  dt2 %>% 
  filter(sub == "h1",
         outcome %in% dths_names) %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))

dts_h3 <- 
  dt2 %>% 
  filter(sub == "h3") %>% 
  summarise(value = n(), 
            .by = c(year, sex, age))




# 
flu0 %>% 
  filter(sub %in% c("h1", "h3")) %>% 
  ggplot()+
  geom_line(aes(cohort, value, col = sub))+
  theme_bw()

# flu circulation by sex, age, and subtype ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flu2 <- 
  flu %>% 
  # filter(sub %in% c("h1", "h3", "b")) %>% 
  # mutate(age2 = interval(ymd(date_bth), ymd(date_flu)) %>% as.numeric('years'),
  #        age2 = round(age2),
  #        age = ifelse(is.na(age), age2, age)) %>% 
  # select(-age2) %>% 
  summarise(value = n(), .by = c(year, sex, age, sub, hosp, outcome)) %>% 
  drop_na(age)

ili <- 
  flu2 %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "cases")

css <- 
  flu2 %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "cases")

hsp <- 
  flu2 %>% 
  filter(hosp == 1) %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "hosps")

dts <- 
  flu2 %>% 
  filter(outcome == "death_flu") %>% 
  summarise(value = sum(value), .by = c(year, sex, age, sub)) %>% 
  mutate(outcome = "deaths")

flu3 <- 
  bind_rows(dts, css) %>% 
  group_by(year, age, outcome, sub) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(dts, css) %>% 
  mutate(cohort = year - age) %>% 
  left_join(pop2) %>% 
  drop_na(age) %>% 
  mutate(mx = 1e5 * value / pop)

flu3 %>% 
  filter(sex == "t") %>% 
  group_by(sub, outcome) %>% 
  summarise(sum(value))

write_rds(flu3, "data_inter/flu_data_brazil_exposures_2009_2019.rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flu3 <- read_rds("data_inter/flu_data_brazil_exposures_2009_2019.rds")

# H1 ====
# ~~~~~~~
h1 <- 
  flu3 %>% 
  filter(sub == "h1",
         sex == "t",
         outcome == "deaths") %>% 
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
  flu3 %>% 
  filter(sub == "h3",
         sex == "t",
         outcome == "death") %>% 
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
  flu3 %>% 
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

