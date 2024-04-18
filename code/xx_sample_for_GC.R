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


chunk <- 
  all %>% 
  filter(year == "2016", 
         type == "h1",
         # cause == "pi",
         sex == "t") %>% 
  mutate(age = ifelse(age >= 90, 90, age)) %>% 
  reframe(hsp = sum(hsp),
          dts = sum(dts),
          exposure = sum(pop),
          .by = c(age))

write_rds(chunk, "data_inter/sample_h1_2016.rds")


chunk %>% 
  ggplot()+
  geom_line(aes(age, dts, col = year))+
  scale_y_log10()
