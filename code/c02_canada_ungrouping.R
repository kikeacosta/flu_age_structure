rm(list=ls())
library(tidyverse)
library(ungroup)
options(scipen=999)

ca <- read_csv("data_input/canada_2023w05.csv")
au <- read_rds("data_inter/subtypes_age_australia_std.rds")
pop_ca <- read_csv("data_input/Canada_exposure_hmd.csv")
pop_au <- read_csv("data_input/Australia_exposure_hmd.csv")


# population data
pop_au2 <- 
  pop_au %>% 
  filter(Year >= 2000) %>% 
  select(year = Year, age = Age, pop = Total) %>% 
  mutate(age = ifelse(age == "110+", 110, age %>% as.integer()),
         age  = ifelse(age > 100, 100, age)) %>% 
  group_by(year, age) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()



# processing Canadian data and imputing h1 and h3 flu cases
ca2 <- 
  ca %>% 
  select(age = 1,
         prop_a = 2,
         cases = 4,
         prop_h1 = 5,
         prop_h3 = 6,
         cases_a_sub = 7) %>% 
  mutate(prop_h1 = prop_h1 %>% str_replace("%", "") %>% as.double() / 100,
         prop_h3 = prop_h3 %>% str_replace("%", "") %>% as.double() / 100,
         prop_a = prop_a %>% str_replace("%", "") %>% as.double() / 100,
         cases_a = cases * prop_a,
         cases_h1_sub = cases_a_sub * prop_h1,
         cases_h3_sub = cases_a_sub * prop_h3,
         h1 = cases_a*cases_h1_sub/(cases_h1_sub+cases_h3_sub),
         h3 = cases_a*cases_h3_sub/(cases_h1_sub+cases_h3_sub)) %>% 
  select(age, h1, h3) %>% 
  mutate(age = str_sub(age, 1, 2) %>% str_trim() %>% as.integer())

unique(ca2$age)


# Australia data
# ~~~~~~~~~~~~~~
# converting Australia data to Canada ages to see how well ungrouping works
au2 <- 
  au %>% 
  mutate(age = case_when(age %in% 0:4 ~ 0,
                         age %in% 5:19 ~ 5,
                         age %in% 20:44 ~ 20,
                         age %in% 45:64 ~ 45,
                         age >= 65 ~ 65)) %>% 
  group_by(year, subtype, age) %>% 
  summarise(cases = sum(cases)) %>% 
  ungroup()

chunk <- au2 %>% 
  filter(year == 2009,
         subtype == "h1")

# ungrouping ages in single-year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ungroup_this <- 
  function(chunk,
           lambda = 100){
    
    ages <- chunk$age
    # increasing numbers to avoid inconsistent outcomes (decreasing afterwards)
    cass <- chunk$cases*1e5
    # closing at age 100+
    nlast <- 101 - max(ages)
    yr <- unique(chunk$year)
    pops <- 
      pop_au2 %>% 
      filter(year == yr) %>% 
      pull(pop)
    
    V1 <- pclm(x = ages,
               y = cass,
               nlast = nlast)$fitted
    V2 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               control = list(lambda = lambda, deg = 3))$fitted
    V3 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               offset = pops)$fitted * pops
    V4 <- pclm(x = ages,
               y = cass,
               nlast = nlast,
               offset = pops,
               control = list(lambda = lambda, deg = 3))$fitted * pops
    
    out <- 
      tibble(age = 0:100, 
             cases1 = V1/1e5,
             cases2 = V2/1e5,
             cases3 = V3/1e5,
             cases4 = V4/1e5)
    return(out)
  }

au3 <- 
  au2 %>% 
  arrange(year, subtype, age) %>% 
  group_by(year, subtype) %>% 
  do(ungroup_this(chunk = .data, lambda = 1e10)) %>% 
  ungroup()

au4 <- 
  au3 %>% 
  mutate(age = age-age%%5,
         age = ifelse(age >= 85, 85, age)) %>% 
  gather(starts_with("cases"), key = method, value = cases) %>% 
  group_by(year, age, subtype, method) %>% 
  summarise(cases = sum(cases)) %>% 
  ungroup()

au5 <- 
  au %>% 
  mutate(method = "original") %>% 
  bind_rows(au4)

au5 %>% 
  # filter(year == 2015) %>% 
  ggplot()+
  geom_point(aes(age, cases, col = method))+
  facet_wrap(year~subtype, scales = "free_y", ncol = 2)+
  theme_bw()
ggsave("figures/ungrouping_test_australia_using_canada_ages.pdf",
       w = 5,
       h = 30)


au5 %>% 
  filter(method %in% c("original", "cases4")) %>%
  ggplot()+
  geom_point(aes(age, cases, col = method))+
  facet_wrap(year~subtype, scales = "free_y", ncol = 2)+
  theme_bw()
ggsave("figures/ungrouping_test_australia_using_canada_ages2.pdf",
       w = 5,
       h = 30)


au5 %>% 
  filter(method %in% c("original", "cases4"),
         year == 2009,
         subtype == "h3") %>%
  ggplot()+
  geom_point(aes(age, cases, col = method))+
  facet_wrap(year~subtype, scales = "free_y", ncol = 2)+
  theme_bw()
ggsave("figures/example_ungrouping_working.png",
       w = 4,
       h = 3)


# comparing average ages
test <- 
  au5 %>% 
  mutate(age2 = case_when(age %in% 0:4 ~ 0,
                          age %in% 5:19 ~ 5,
                          age %in% 20:44 ~ 20,
                          age %in% 45:64 ~ 45,
                          age >= 65 ~ 65)) %>% 
  group_by(year, method, subtype, age2) %>% 
  summarise(av_age = sum((age+2.5)*cases)/sum(cases)) %>% 
  ungroup()
