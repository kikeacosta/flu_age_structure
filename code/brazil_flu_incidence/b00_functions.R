rm(list=ls())
library(wpp2022)
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
library(mgcv)

options(scipen=999)

# smoothing rates over ages
smooth_age <- function(chunk){
  try(
    model <- 
    gam(value ~ 
          s(age, bs = 'ps', m = c(2,2), k = 30) + 
          offset(log(exposure)), 
        # weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"), gamma = 1e-10)
  )
  
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
      mutate(val_smt = res$fit,
             val_smt_lc = val_smt - 1.96 * res$se.fit,
             val_smt_uc = val_smt + 1.96 * res$se.fit)
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(val_smt = NA,
             val_smt_lc = NA,
             val_smt_uc = NA)
  }
  return(chunk2)
}
