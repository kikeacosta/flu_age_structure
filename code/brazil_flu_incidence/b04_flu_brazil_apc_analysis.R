rm(list=ls())
{
  libs <- c("Epi", "MASS", "writexl", "tidyverse")
  for (i in libs){
    library(i,character.only = TRUE)
  }
  select <- dplyr::select
}
source("code/00_functions.R")


# loading data ====
# ~~~~~~~~~~~~~~~~~

# weekly positive influenza cases by subtype, sex, and age
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cs <- 
  read_rds("data_inter/flu_data_brazil_exposures_2009_2019.rds")

# transforming it to wide format, so we can chose numerator and denominator
# (death rates, incidence, case fatality rates)

db_in <- 
  cs %>% 
  mutate(value = round(value) + 1) %>% 
  select(-mx) %>% 
  spread(outcome, value) %>% 
  replace_na(list(cases = 0,
                  deaths = 0))

# sex <- "t"
# sub = "h1"
# num <- "cases"
# den <- "population"

{ # selecting group sizes, and age and period limits
  # size of categories
  gr <- 1
  # ages
  amin <- 0;
  amax <- 85;
  # periods (maximum 2014, because the first group includes seasonal years 1959-60 & 1960-61, so
  # the group includes the seasonal periods 2013-14 & 2014-15, thus, seasonal period 2015-16 is excluded)
  pmin <- 2009; 
  pmax <- 2019; 

  # distribution: negative binomial
  dist <- "nb"
  # drift extraction (consult Crastensen material), here, naive (n)
  extr <- "n" 
}

parm="APC"; mod="factor";
dr.extr="n";
amn=amin; amx=amax;
pmn=pmin; pmx=pmax;
dis=dist

# virus subtype
sb <- sub
# numerator
nm <- "cases"
# denominator
dn <- "population"
# sex
sx <- "t"

# The function for fitting APC models (apc_acp) can be found in the "00_functions.R" script

# cohort effects for incidence, mortality and CFR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  # on incidence
  # ~~~~~~~~~~~~
  apc_h1_incid <- 
    apc_acp(sx = "t", sb = "h1", 
            nm = "cases", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  # it has two components: "estimates" and "plot"
  
  # plotting cohort effects "plot"
  apc_h1_incid$plot
  
  # extracting other coefficients, it can be: "Age", "Period", or "Cohort"
  h1_incid_effects <- 
    apc_h1_incid$estimates %>% 
    filter(dimension == "Age")
  
  # plotting them
  h1_incid_effects %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    scale_y_log10()+
    theme_bw()
  
  
  # on mortality by h1n1
  # ~~~~~~~~~~~~~~~~~~~~
  apc_h1_death <- 
    apc_acp(sx = "t", sb = "h1", 
            nm = "deaths", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  apc_h1_death$plot
  
  # case fatality rate (CFR)
  # ~~~~~~~~~~~~~~~~~~~~~~~~
  apc_h1_cfr <- 
    apc_acp(sx = "t", sb = "h1", 
            nm = "deaths", dn = "cases",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  
  effects <- 
    apc_h1_cfr$estimates %>% 
    filter(dimension == "Age")
  
  test %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    scale_x_continuous(breaks = seq(0, 100, 10))
  
  
  apc_h1_cfr$plot
  
  # incidence
  apc_h3_incid <- 
    apc_acp(sx = "t", sb = "h3", 
            nm = "cases", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  # mortality
  apc_h3_death <- 
    apc_acp(sx = "t", sb = "h3", 
            nm = "deaths", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  # cfr
  apc_h3_cfr <- 
    apc_acp(sx = "t", sb = "h3", 
            nm = "deaths", dn = "cases",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)

  
  
} 
  
# ratios H1/H3 for cohort effects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
apc_incid <- 
  apc_h1_incid[[1]] %>% 
  select(Years, rr_h1 = value, dimension) %>% 
  left_join(apc_h3_incid[[1]] %>% 
              select(Years, rr_h3 = value, dimension)) %>% 
  mutate(rh1h3 = rr_h1 / rr_h3)

apc_incid %>% 
  filter(dimension == "Cohort",
         Years %in% 1930:2016) %>% 
  ggplot()+
  geom_line(aes(Years, rh1h3))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
  annotate("text", 
           y = 1, 
           x = c(1947, 1957, 1968, 1978, 2009), 
           label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
           angle = 90, hjust = 1, vjust = 0, size = 3)+
  scale_x_continuous(breaks = seq(1910, 2010, 10))+
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5))+
  labs(x = "Cohort", title = "ratio H1N1 / H3N2")+
  theme_bw()


  apc_h1_death[[1]] %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    scale_y_log10()+
    facet_grid(~dimension, scales = "free", space = "free")+
    geom_hline(yintercept = 1, linetype = "dashed")+
    theme_bw()


  apc_h1_cfr$estimates %>% 
    filter(dimension == "Age") %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_x_continuous(breaks = seq(0, 100, 10))+
    scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5))+
    labs(x = "Age", title = "H1N1 CFR")+
    theme_bw()
  

detach(package:MASS)
detach(package:Epi)
