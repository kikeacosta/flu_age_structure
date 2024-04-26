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
source("code/00_functions.R")

# loading influenza data (cases, hospitalizations, deaths) in Brazil 
dt_in <- 
  readRDS("data_inter/sample_dts_hsp_bra.rds")

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
    30, 
    # age interval
    ag1 = 0, 
    ag2 = 90, 
    # year
    2009, 
    # type of flu
    "h1")

# plot of splines and slopes
t09$p_psplines
# plot of splines and slopes
t09$p_sizer

# estimates of p-splines (type1 indicates whether are Log mx estimates of slopes)
t09$dt_ests
# data for SiZer
t09$dt_sizer


# testing different knots
for(kn in c(5, 10, 18, 19, 20, 30, 50, 75, 100, 200, 500)){
  do_gc_magic(dt_in, kn, ag1 = 0, ag2 = 90, 2016, "h1")
}

## END
