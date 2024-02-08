rm(list=ls())
source("code/brazil_flu_incidence/b00_functions.R")

pop <- 
  read_rds("data_inter/brazil_exposures_1999_2022_ages_0_110.rds")

dt <- 
  read_rds("data_inter/flu_data_brazil_2009_2019_v3.rds") %>% 
  mutate(age = year - cohort)

{
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
}

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


# # deaths
# dt2 %>% 
  



# descriptive plot of data
# ~~~~~~~~~~~~~~~~~~~~~~~~
cols <- 
  c("h1" = "#d00000", 
    "h3" = "#3a86ff")

css_epi <- 
  dt2 %>% 
  filter(sub %in% c("h1", "h3"),
         sex == "t",
         flu == 1) %>% 
  summarise(value = sum(value), .by = c(year, sub)) %>% 
  mutate(type = "Cases")

dts_epi <- 
  dt2 %>% 
  filter(sub %in% c("h1", "h3"),
         sex == "t",
         flu == 1) %>% 
  mutate(value = value*dth) %>% 
  summarise(value = sum(value), .by = c(year, sub)) %>% 
  mutate(type = "Deaths")

dts_epi %>% summarise(value = sum(value), .by = sub)

props <- 
  dts_epi %>% 
  spread(sub, value) %>% 
  mutate(t = h1+h3,
         h1 = h1/t, 
         h3 = h3/t) %>% 
  select(-t) %>% 
  gather(h1, h3, key = sub, value = value) %>% 
  mutate(type = "Share")

css_dts_epi <- 
  bind_rows(css_epi, dts_epi, props)

css_dts_epi %>% 
  ggplot(aes(fill=sub, y=value, x=year))+
  facet_wrap(~type, scales = "free_x")+
  geom_bar(position="stack", stat="identity")+
  scale_x_continuous(breaks = 2009:2019)+
  scale_y_continuous()+
  scale_fill_manual(values = cols)+
  theme_bw()+
  theme(strip.background = element_blank(),
        legend.position = "left")+
  coord_flip(expand = 0)+
  labs(y = "Flu cases", x = "Year", fill = "Flu\nsubtype")

# ggsave("figures/brazil/bcn_slides/surv_data_cases_deaths.png",
#        w = 8, h = 3)

css_epi %>% 
  ggplot(aes(fill=sub, y=value, x=year))+
  geom_bar(position="fill", stat="identity")+
  scale_x_continuous(breaks = 2009:2019)+
  # scale_y_continuous(limits = c(0, 15000))+
  scale_fill_manual(values = cols)+
  theme_bw()+
  theme(legend.position = "left")+
  coord_flip(expand = 0)+
  labs(y = "Flu deaths", x = "Year", fill = "Flu\nsubtype")

# ggsave("figures/brazil/bcn_slides/surv_data_cases_prop.png",
#        w = 6, h = 3)



# grouping data 
{
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
  filter(
    sex == "t",
    # flu != 1,
    # sub == "z",
    !(sub %in% c("h1", "ab"))
    ) %>% 
  filter() %>% 
  # mutate(dth = ifelse(outcome == "death_flu", 1, 0)) %>% 
  summarise(css = sum(value),
            hsp = sum(hosp*value),
            dts = sum(dth*value),
            .by = c(year, sex, age)) %>% 
  mutate(type = "nonflu")
}
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
{
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
}

# plots ====
# ~~~~~~~~~~
# Incidence

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

# ggsave("figures/brazil/bcn_slides/surv_ix_h1_age.png",
#        w = 8, h = 4)


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
# ggsave("figures/brazil/bcn_slides/surv_ix_h1_coh.png",
#        w = 8, h = 4)

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
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed")

# ggsave("figures/brazil/bcn_slides/surv_ix_h1_coh_lines.png",
#        w = 8, h = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mortality ====
# ~~~~~~~~~~~~~~
# death rates
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
# ggsave("figures/brazil/bcn_slides/surv_mx_h1_age.png",
#        w = 8, h = 4)


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
# ggsave("figures/brazil/bcn_slides/surv_mx_h1_coh.png",
#        w = 8, h = 4)


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
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed", 
             col = c("#7b2cbf", "#3a86ff"))+
  geom_text(aes(x = 1969, y = 0.02),
            label = "1968 H3N2 pandemic", angle = 90,
            size = 3, hjust = 0, col = "#3a86ff")+
  geom_text(aes(x = 1958, y = 0.02), 
            label = "1957 H2N2 pandemic", angle = 90, 
            size = 3, hjust = 0, col = "#7b2cbf")

# ggsave("figures/brazil/bcn_slides/surv_mx_h1_coh_lines.png",
#        w = 8, h = 4)


mx %>% 
  mutate(cohort = year - age,
         value_r = 1e5*value / exposure,
         val_smt_r = 1e5*val_smt/exposure) %>% 
  filter(sex == "t",
         type %in% c("h1", "nonflu"),
         year %in% c(2009),
         cohort %in% 1925:2010) %>%
  ggplot()+
  geom_line(aes(cohort, val_smt_r, linetype = type),
            col = "#6a040f")+
  # geom_point(aes(cohort, value_r, col = factor(year), linetype = type), 
  #            alpha = 0.4, size = 0.5)+
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  # scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(y = "death rates (/100K)")+
  theme_bw()+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed", 
             col = c("#7b2cbf", "#3a86ff"))+
  geom_text(aes(x = 1969, y = 0.02),
            label = "1968 H3N2 pandemic", angle = 90,
            size = 3, hjust = 0, col = "#3a86ff")+
  geom_text(aes(x = 1958, y = 0.02), 
            label = "1957 H2N2 pandemic", angle = 90, 
            size = 3, hjust = 0, col = "#7b2cbf")

# ggsave("figures/brazil/bcn_slides/surv_mx_h1_ref_coh_lines.png",
#        w = 8, h = 4)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# relative risks ====
# ~~~~~~~~~~~~~~~~~~~
# H1 relative to all not identified as flu
unique(mx$type)

rrh1 <- 
  mx %>% 
  filter(type == "h1") %>% 
  select(year, age, css_h1 = css, dts_h1 = value, dts_h1_smt = val_smt) %>% 
  left_join(mx %>% 
                filter(type == "nonflu") %>% 
                select(year, age, css_nf = css, dts_nf = value, dts_nf_smt = val_smt)) %>% 
  mutate(rr_h1 = dts_h1/dts_nf,
         rr_h1_smt = dts_h1_smt/dts_nf_smt,
         cohort = year - age)


rrh1 %>% 
  filter(
    # sex == "t",
    # type == "h1",
    year %in% c(2009),
    cohort %in% 1925:2010
  ) %>%
  ggplot()+
  geom_line(aes(cohort, rr_h1_smt, col = factor(year)))+
  geom_point(aes(cohort, rr_h1, col = factor(year)), alpha = 0.4, size = 0.5)+
  
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(col = "Year", y = "relative risks")+
  theme_bw()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed", 
             col = c("#7b2cbf", "#3a86ff"))+
  geom_text(aes(x = 1969, y = 0.02),
            label = "1968 H3N2 pandemic", angle = 90,
            size = 3, hjust = 0, col = "#3a86ff")+
  geom_text(aes(x = 1958, y = 0.02), 
            label = "1957 H2N2 pandemic", angle = 90, 
            size = 3, hjust = 0, col = "#7b2cbf")


ggsave("figures/brazil/bcn_slides/surv_rr_h1_pan_coh_lines.png",
       w = 8, h = 4)

rrh1 %>% 
  filter(
    # sex == "t",
    # type == "h1",
    year %in% c(2009, 2013, 2016),
    cohort %in% 1925:2010
  ) %>%
  ggplot()+
  geom_line(aes(cohort, rr_h1_smt, col = factor(year)))+
  geom_point(aes(cohort, rr_h1, col = factor(year)), alpha = 0.4, size = 0.5)+
  
  scale_y_log10()+
  scale_x_reverse(breaks = seq(1900, 2010, 10))+
  scale_color_manual(values = cols)+
  # facet_wrap(~year, scales = "free_y")+
  labs(col = "Year", y = "relative risks")+
  theme_bw()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = c(1957, 1968), linetype = "dashed", 
             col = c("#7b2cbf", "#3a86ff"))+
  geom_text(aes(x = 1969, y = 0.02),
            label = "1968 H3N2 pandemic", angle = 90,
            size = 3, hjust = 0, col = "#3a86ff")+
  geom_text(aes(x = 1958, y = 0.02), 
            label = "1957 H2N2 pandemic", angle = 90, 
            size = 3, hjust = 0, col = "#7b2cbf")


ggsave("figures/brazil/bcn_slides/surv_rr_h1_coh_lines.png",
       w = 8, h = 4)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

