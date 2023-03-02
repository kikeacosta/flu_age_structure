rm(list=ls())
source("R/00_functions.R")

# loading data ====
# ~~~~~~~~~~~~~~~~~

# weekly positive influenza cases by subtype, sex, and age
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
css0815 <- 
  read_xlsx("data_input/national-notifiable-diseases-surveillance-system-nndss-public-dataset-influenza-laboratory-confirmed-dataset.xlsx",
            skip = 4,
            sheet = 1)

css1618 <- 
  read_xlsx("data_input/national-notifiable-diseases-surveillance-system-nndss-public-dataset-influenza-laboratory-confirmed-dataset.xlsx",
            skip = 4,
            sheet = 2)

css1921 <- 
  read_xlsx("data_input/national-notifiable-diseases-surveillance-system-nndss-public-dataset-influenza-laboratory-confirmed-dataset.xlsx",
            skip = 4,
            sheet = 3)


# preparing cases ====
# ~~~~~~~~~~~~~~~~~~~~

css <- 
  bind_rows(css0815, 
            css1618, 
            css1921)

unique(css$State)

css2 <- 
  css %>% 
  rename(date = 1,
         state = 2,
         age = 3, 
         sex = 4,
         ind = 5,
         sub = 6) %>% 
  mutate(date = ymd(date))

unique(css2$state)
unique(css2$age)
unique(css2$sex)
unique(css2$sub)

css3 <- 
  css2 %>% 
  arrange(date) %>% 
  group_by(date, sex, age, sub) %>% 
  summarise(css = n()) %>% 
  ungroup() %>% 
  spread(sub, css, fill = 0) %>% 
  rename(ab = 4,
         h1pre = 5,
         h1 = 6,
         h3 = 7,
         a_uns = 8, 
         b = 9, 
         c = 10, 
         unt = 11) %>% 
  gather(-date, -sex, -age, key = sub, value = css) 

# totals
# ~~~~~~

css3 %>% 
  group_by(sub) %>% 
  summarise(css = sum(css)) 

css3 %>% 
  group_by(date) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_line(aes(date, css))

css3 %>% 
  group_by(date, age) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_line(aes(date, css))+
  facet_wrap(~age, scales = "free_y")

css3 %>% 
  group_by(date, sex) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_line(aes(date, css, col = sex))

css3 %>% 
  group_by(date, sub) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  ggplot()+
  geom_line(aes(date, css, col = sub))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subtype imputation ====
# ~~~~~~~~~~~~~~~~~~~~~~~

tots <- 
  css3 %>% 
  group_by(date, sex, age) %>% 
  summarise(tot = sum(css)) %>% 
  ungroup()

css4 <- 
  css3 %>% 
  filter(sub != "unt") %>% 
  left_join(tots) %>% 
  group_by(date, sex, age) %>% 
  mutate(css = ifelse(sum(css) * tot == 0, 
                      0,
                      tot * (css/sum(css)))) %>% 
  ungroup() %>% 
  select(-tot)

unique(css4$sub)

as <- c("a_uns", "h1", "h1pre", "h3")

css_a <- 
  css4 %>% 
  filter(sub %in% as)

a_tot <- 
  css_a %>% 
  group_by(date, sex, age) %>% 
  summarise(a_tot = sum(css)) %>% 
  ungroup()

css_a2 <- 
  css_a %>% 
  filter(sub != "a_uns") %>% 
  left_join(a_tot) %>% 
  group_by(date, sex, age) %>% 
  mutate(css = ifelse(sum(css) * a_tot == 0, 
                      0,
                      a_tot * (css/sum(css)))) %>% 
  ungroup() %>% 
  select(-a_tot)

css4 %>% 
  group_by(sub) %>% 
  summarise(css = sum(css))

css5 <- 
  css4 %>% 
  filter(!sub %in% c(as, "ab", "c")) %>% 
  bind_rows(css_a2) %>% 
  arrange(date, sex, age) %>% 
  mutate(sex = str_sub(sex, 1, 1) %>% str_to_lower())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# age imputation ====
# ~~~~~~~~~~~~~~~~~~~~~~~

tot_age <- 
  css5 %>% 
  group_by(date, sex, sub) %>% 
  summarise(tot_age = sum(css)) %>% 
  ungroup()

unique(css5$age)


css6 <- 
  css5 %>% 
  filter(age != "Unknown") %>% 
  left_join(tot_age) %>% 
  group_by(date, sex, sub) %>% 
  mutate(css = ifelse(sum(css) * tot_age == 0, 
                      0,
                      tot_age * (css/sum(css)))) %>% 
  ungroup() %>% 
  select(-tot_age) 

# sex imputation ====
# ~~~~~~~~~~~~~~~~~~~~~~~

tot_sex <- 
  css6 %>% 
  group_by(date, age, sub) %>% 
  summarise(tot_sex = sum(css)) %>% 
  ungroup()

unique(css6$sex)

css7 <- 
  css6 %>% 
  filter(sex != "u") %>% 
  left_join(tot_sex) %>% 
  group_by(date, sex, sub) %>% 
  mutate(css = ifelse(sum(css) * tot_sex == 0, 
                      0,
                      tot_sex * (css/sum(css)))) %>% 
  ungroup() %>% 
  select(-tot_sex) 


# adding total sex ====
# ~~~~~~~~~~~~~~~~~~~~~
css8 <- 
  css7 %>% 
  group_by(date, age, sub) %>% 
  summarise(css = sum(css)) %>% 
  ungroup() %>% 
  mutate(sex = "t") %>% 
  bind_rows(css7) 

write_rds(css8, "data_inter/weekly_flu_subtypes_australia.rds")



