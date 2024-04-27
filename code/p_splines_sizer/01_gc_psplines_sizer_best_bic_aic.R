## R-code for reproducing SiZer maps using P-splines
## SiZer map consists in plotting sign and significance of the derivative
## of the linear predictor over x-axis and different smoothing parameters

## Application death counts (or hospitalizations)
## by single-year age for H1N1-flu-infected people in 2016 in Brazil.
## In this case, data are smooth over age in a Poisson setting
## therefore the derivative over age of the linear predictor
## corresponds to what demographers call Rate-of-Aging 
## (or LAR, Lifetable Aging Rates) basically the derivatives 
## of the log-mortality or the relative derivative of the force of mortality

## Data provided by Kike Acosta
## by Giancarlo Camarda, 2024.03.20

## Adapted and almost destroyed by Kike Acosta, 2024.03.23
## wrapped in a function

rm(list = ls())
# calling functions built from Giancarlo's code
source("code/p_splines_sizer/00_functions_iter.R")

# loading influenza data (cases, hospitalizations, deaths) in Brazil 
dt_in <- 
  readRDS("data_inter/master_brazil_flu_2009_2023.rds")

unique(dt_in$year)
unique(dt_in$type)

# function to estimate p-splines according to penalization levels (lambda)
# it provides the estimates of the p-splines, the best AIC and BIC
# and plot the splines and SiZer plots of the slopes (save them in figures:/)
t09 <-
  do_gc_magic(
    # data
    dt_in,
    # number of knots
    20,
    # age interval
    ag1 = 0,
    ag2 = 90,
    # year
    2018,
    # type of flu
    "h1")

# plot of splines and slopes
t09$aic
# plot of splines and slopes
t09$bic

# # estimates of p-splines (type1 indicates whether are Log mx
# estimates of slopes)
# t09$dt_ests

dt <- 
  dt_in %>% 
  filter(
    # year == 2018,
         age %in% 0:90,
         type == "h1")

# Iterating across all years, and all possible knots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tp <- "h1"

all_aic_ests <- tibble()
all_bic_ests <- tibble()
all_aic_sizers <- tibble()
all_bic_sizers <- tibble()

for(yr in 2009:2023){
  
  sums <- tibble()
  aic_ests <- tibble()
  bic_ests <- tibble()
  aic_sizer <- tibble()
  bic_sizer <- tibble()
  
  dt <- 
    dt_in %>% 
    filter(
      year == yr,
      age %in% 0:90,
      type == tp
    )
  
  for(k in 5:50){
    
    tmp <- 
      get_best_aic_bic(
        # data
        dt,
        # number of knots
        k,
        # age interval
        ag1 = 0,
        ag2 = 90,
        # year
        yr,
        # type of flu
        tp
      )
    
    sums <- 
      sums %>% 
      bind_rows(
        tibble(year = yr,
               knots = k,
               aic = tmp$aic,
               bic = tmp$bic,
        )
      )
  }
  
  sums_bic <- 
    sums %>% 
    arrange(bic) %>% 
    filter(bic < min(bic) + 6) %>% 
    slice_head(n = 3)
  
  sums_aic <- 
    sums %>% 
    arrange(aic) %>% 
    filter(aic < min(aic) + 6) %>% 
    slice_head(n = 3)
  
  best_aic_ks <- 
    sums_aic %>% 
    pull(knots)
  
  best_bic_ks <- 
    sums_bic %>% 
    pull(knots)
  
  # i <- 1
  for(i in 1:3){
    
    tmp_aic <- 
      do_gc_magic_best(
        # data
        dt,
        # number of knots
        best_aic_ks[i],
        # age interval
        ag1 = 0,
        ag2 = 90,
        # year
        yr,
        # type of flu
        "h1"
      )
    
    tmp_bic <- 
      do_gc_magic_best(
        # data
        dt,
        # number of knots
        best_bic_ks[i],
        # age interval
        ag1 = 0,
        ag2 = 90,
        # year
        yr,
        # type of flu
        "h1"
      )
    
    # Best p-splines
    aic_ests <- 
      aic_ests %>% 
      bind_rows(
        tmp_aic$dt_ests %>% 
          filter(opt == "AIC") %>% 
          mutate(stat = "aic",
                 best_pos = i,
                 year = yr)
      )
    
    bic_ests <- 
      bic_ests %>% 
      bind_rows(
        tmp_bic$dt_ests %>% 
          filter(opt == "BIC") %>% 
          mutate(stat = "bic",
                 best_pos = i,
                 year = yr)
      )
    
    # first derivative of splines for SiZer construction
    
    aic_sizer <- 
      aic_sizer %>% 
      bind_rows(
        tmp_aic$dt_sizer %>% 
          filter(opt == "AIC") %>% 
          mutate(stat = "aic",
                 best_pos = i,
                 year = yr)
      )
    
    bic_sizer <- 
      bic_sizer %>% 
      bind_rows(
        tmp_bic$dt_sizer %>% 
          filter(opt == "BIC") %>% 
          mutate(stat = "bic",
                 best_pos = i,
                 year = yr)
      )
    
  }
  
  all_aic_ests <- 
    all_aic_ests %>% 
    bind_rows(aic_ests)
  
  all_bic_ests <- 
    all_bic_ests %>% 
    bind_rows(bic_ests)
  
  all_aic_sizers <- 
    all_aic_sizers %>% 
    bind_rows(aic_sizer)
  
  all_bic_sizers <- 
    all_bic_sizers %>% 
    bind_rows(bic_sizer)
  
}

all_aic_sizers %>%
  ggplot(aes(cohort, best_pos)) + 
  facet_grid(year~., switch = "y")+
  geom_tile(aes(fill = sign))+
  scale_fill_manual(name = NULL,
                    values = c("0"="#DDCC77", "-1"="#88CCEE", "1"="#CC6677"),
                    labels = c("Significant Negative", "Not Significant", "Significant Positive"))+
  geom_vline(xintercept = c(1957, 1968, 1978, 2009), linetype = "dashed",
             col = "grey40")+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  coord_cartesian(xlim = c(1919, 2023),
                  expand = 0)+
  labs(title = "SiZer with best AIC by year",
       y = "Year", x = "Cohort")+
  theme_minimal() +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y  = element_blank(),
        legend.text = element_text(size = 12),
        legend.position="top",
        panel.spacing = unit(0.1,'lines'),
        plot.title = element_text(size=18, hjust = 0),
        plot.subtitle = element_text(size=16, hjust = 0),
        strip.placement = "outside"
  )

ggsave(paste0("figures/sizers_best/", tp, "_sizer_best_aics_2009_2023.png"),
       w = 10, h = 10)

ggsave(paste0("figures/sizers_best/", tp, "_sizer_best_aics_2009_2023.pdf"),
       w = 10, h = 10)

all_bic_sizers %>%
  ggplot(aes(cohort, best_pos)) + 
  facet_grid(year~., switch = "y")+
  geom_tile(aes(fill=sign))+
  scale_fill_manual(name = NULL,
                    values = c("0"="#DDCC77", "-1"="#88CCEE", "1"="#CC6677"),
                    labels = c("Significant Negative", "Not Significant", "Significant Positive"))+
  geom_vline(xintercept = c(1957, 1968, 1978, 2009), linetype = "dashed",
             col = "grey40")+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  coord_cartesian(xlim = c(1919, 2023),
                  expand = 0)+
  labs(title = "SiZer with best BIC by year",
        y = "Year", x = "Cohort")+
  theme_minimal() +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        axis.text.y  = element_blank(),
        legend.text = element_text(size = 12),
        legend.position="top",
        panel.spacing = unit(0.1,'lines'),
        plot.title = element_text(size=18, hjust = 0),
        plot.subtitle = element_text(size=16, hjust = 0),
        strip.placement = "outside"
  )

ggsave(paste0("figures/sizers_best/", tp, "_sizer_best_bics_2009_2023.png"),
       w = 10, h = 10)

ggsave(paste0("figures/sizers_best/", tp, "_sizer_best_bics_2009_2023.pdf"),
       w = 10, h = 10)

