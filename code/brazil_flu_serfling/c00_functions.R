rm(list=ls())
library(tidyverse)
library(lubridate)
library(ungroup)
library(readxl)
library(viridis)
library(mgcv)
options(scipen=999)
set.seed(2019) 

# smoothing rates over ages
smooth_age <- function(chunk){
  model <- 
    gam(flu ~ 
          s(age, bs = 'ps', m = c(2,2), k = 30) + 
          offset(log(exposure)), 
        # weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"), gamma = 1e-10)
  
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
      mutate(flu_smt = res$fit,
             flu_smt_lc = bsn - 1.96 * res$se.fit,
             flu_smt_uc = bsn + 1.96 * res$se.fit)
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(flu_smt = NA,
             flu_smt_lc = NA,
             flu_smt_uc = NA)
  }
  return(chunk2)
}


# monthly baselines
# ~~~~~~~~~~~~~~~~~
est_mth_epis <- function(chunk){
  
  model <- 
    gam(dts ~ 
          s(t, bs = 'ps', m = c(2,2)) + 
          sn2 + cs2 +
          sn4 + cs4 +
          sn8 + cs8 +
          sn10 + cs10 +
          offset(log(exposure)), 
        weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"))
  
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
      mutate(bsn = res$fit,
             bsn_lc = bsn - 1.96 * res$se.fit,
             bsn_uc = bsn + 1.96 * res$se.fit)
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(bsn = NA,
             bsn_lc = NA,
             bsn_uc = NA)
  }
  
  return(chunk2)
}

est_mth_baseline_pi <- function(chunk){
  
  model <- 
    gam(dts ~ 
          s(t, bs = 'ps', m = c(2,2)) + 
          sn2 + cs2 +
          sn4 + cs4 +
          sn8 + cs8 +
          sn10 + cs10 +
          offset(log(exposure)), 
        weights = w,
        data = chunk, 
        family = quasipoisson(link = "log"))
  
  test <- 
    try(
      res <- 
        predict(model, 
                newdata = chunk %>% mutate(w = 1),
                type = "response", 
                se.fit = TRUE)
    )
  
  try(
    chunk2 <- 
      chunk %>% 
      mutate(bsn = res$fit,
             bsn_lc = bsn - 1.96 * res$se.fit,
             bsn_uc = bsn + 1.96 * res$se.fit) %>% 
      left_join(simul_mth_intvals(model, 
                                  model_type = "gam", 
                                  db = chunk, 
                                  nsim = 1000,
                                  p = 0.95),
                by = "t")
  )
  
  if(class(test) == "try-error"){
    chunk2 <- 
      chunk %>% 
      mutate(bsn = NA,
             bsn_lp = NA,
             bsn_up = NA,
             bsn_lc = NA,
             bsn_uc = NA)
  }
  
  return(chunk2)
}

simul_mth_intvals <- 
  function(
    # fitted model 
    model, 
    # either GLM or GAM (needed for model matrix extraction step)
    model_type, 
    # prediction data
    db, 
    # number of iterations
    nsim, 
    # prediction intervals' uncertainty level (between 0 and 1)
    p
  ){
    
    # defining upper and lower prediction quantiles
    lp <- (1 - p) / 2
    up <- 1 - lp
    
    # matrix model extraction
    if(model_type == "glm"){
      X_prd <- model.matrix(model, data = db, na.action = na.pass)
    }
    if(model_type == "gam"){
      X_prd <- predict(model, newdata = db, type = 'lpmatrix')
    }
    
    # estimated coefficients
    beta <- coef(model)
    
    # offsets extracted directly from the prediction data
    offset_prd <- matrix(log(db$exposure))
    # model.offset(x)
    
    # extracting variance covariance matrix
    beta_sim <- MASS::mvrnorm(nsim, 
                              coef(model), 
                              suppressWarnings(vcov(model)))
    
    # simulation process
    Ey_sim <- apply(beta_sim, 1, FUN = function (b) exp(X_prd %*% b + offset_prd))
    
    y_sim <- apply(Ey_sim, 2, FUN = function (Ey) {
      y <- mu <- Ey
      # NA's can't be passed to the simulation functions, so keep them out
      idx_na <- is.na(mu) 
      mu_ <- mu[!idx_na] 
      N <- length(mu_)
      phi <- suppressWarnings(summary(model)$dispersion)
      # in case of under-dispersion, sample from Poisson
      if (phi < 1) { phi = 1 }
      y[!idx_na] <- rnbinom(n = N, mu = mu_, size = mu_/(phi-1))      
      return(y)
    })
    
    # from wide to tidy format
    ints_simul <- 
      db %>% 
      select(t)
    
    colnames_y_sim <- paste0('deaths_sim', 1:nsim)
    
    ints_simul[,colnames_y_sim] <- y_sim
    
    # prediction intervals output
    ints_simul <-
      ints_simul %>%
      pivot_longer(cols = starts_with('deaths_sim'),
                   names_to = 'sim_id', values_to = 'deaths_sim') %>%
      group_by(t) %>%
      summarise(
        bsn_lp = quantile(deaths_sim, lp, na.rm = TRUE),
        bsn_up = quantile(deaths_sim, up, na.rm = TRUE), 
        .groups = 'drop'
      ) 
    
    return(ints_simul)
  }
