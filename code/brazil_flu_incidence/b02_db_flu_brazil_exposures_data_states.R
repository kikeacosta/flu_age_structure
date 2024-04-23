rm(list=ls())
source("code/brazil_flu_incidence/b00_functions.R")

path <- "data_input/brazil/projecoes_2018_populacao_idade_simples_2010_2060_20201209.xls"
shts <- excel_sheets(path = path)
# s <- 13
dt1 <- tibble()
for(s in 1:33){
  
  cat(paste0(shts[s], "\n"))
  
  tmp <- read_xls("data_input/brazil/projecoes_2018_populacao_idade_simples_2010_2060_20201209.xls",
                  sheet = s)
  
  tmp2 <- 
    tmp %>% 
    select(1:52) %>% 
    rename(age = 1) %>% 
    mutate(sex = case_when(str_detect(age, "HOMENS - IDADES SIMPLES") ~ "m",
                           str_detect(age, "MULHERES - IDADES SIMPLES") ~ "f",
                           str_detect(age, "TOTAL - IDADES SIMPLES") ~ "t",
                           TRUE ~ NA_character_)) %>% 
    fill(sex) %>% 
    select(age, sex, everything()) %>% 
    drop_na() %>% 
    rename_at(vars(paste0("...", 2:52)), ~ c(2010:2060) %>% as.character()) %>% 
    filter(age != "IDADE") %>% 
    gather(-age, -sex, key = year, value = pop) %>% 
    mutate(iso2 = str_trim(shts[s]))

  dt1 <- 
    dt1 %>% 
    bind_rows(tmp2)
  
}

unique(dt1$iso2)

stt_codes <- 
  read_csv("data_input/brazil/states_codes.csv", 
           locale = readr::locale(encoding = "latin1")) %>% 
  select(cod, iso2)

dt2 <- 
  dt1 %>% 
  left_join(stt_codes) %>% 
  drop_na(cod) %>% 
  filter(age != "TOTAL") %>% 
  mutate(age = ifelse(age == "90+", 90, age %>% as.double()),
         year = year %>% as.double())

dt3 <- 
  dt2 %>% 
  filter(year == 2010) %>% 
  mutate(year = 2009) %>% 
  bind_rows(dt2) %>% 
  filter(year %in% 2009:2023) %>% 
  select(-cod)

dt3 %>% 
  filter(iso2 == "SP",
         age == 20,
         sex == "t") %>% 
  ggplot()+
  geom_line(aes(year, pop))

write_rds(dt3, "data_inter/brazil_exposures_2009_2023_ages_0_90.rds")

dt3 <- read_rds("data_inter/brazil_exposures_2009_2023_ages_0_90.rds")

ungroup_open_age <- 
  function(chunk){
    
    ages <- chunk$age
    pops <- chunk$pop
    
    nlast <- 101 - max(ages)
    yr <- unique(chunk$year)
    sx <- unique(chunk$sex)
    
    V1 <- pclm(x = ages,
               y = pops,
               nlast = nlast)$fitted
    
    new_pop <- tibble(year = yr,
                      sex = sx,
                      age = 0:100, 
                      pop = V1)
    
    chunk %>% 
      filter(age < 90) %>% 
      bind_rows(new_pop %>% filter(age >= 90)) %>% 
      fill(iso2)
  }

pop2 <- 
  dt3 %>% 
  # filter(iso2 == "SP") %>% 
  group_by(iso2, year, sex) %>% 
  do(ungroup_open_age(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(pop = round(pop))

write_rds(pop2, 
          "data_inter/brazil_exposures_2009_2023_ages_0_100.rds",
          compress = "xz")

