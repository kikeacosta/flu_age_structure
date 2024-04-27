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
library(tidyverse)
library(scales)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("timriffe/MortalitySmooth")
library(MortalitySmooth)

## R-studio to get the same dir as for the .R file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(scipen=999)
set.seed(2019) 

## simple function for "empirical" derivatives
## by GC ~2004 from Stewart's Calculus, sec. 2.7
emp.der <- function(x, y){
  m <- length(x)
  a <- c(diff(y), 1)
  b <- c(diff(x), 1)
  ab <- a/b
  wei <- c(rep(1, m-1), 0)
  
  a1 <- c(1, -1*diff(y))
  b1 <- c(1, -1*diff(x))
  ab1 <- a1/b1
  wei1 <- c(0, rep(1, m-1))
  
  y1emp <- (ab*wei + ab1*wei1)/(wei+wei1)
  return(y1emp)
}

## function to build up B-splines and associated bases for derivatives
BsplineGrad <- function(x, xl, xr, ndx=NULL, deg, knots=NULL){
  if(is.null(knots)){
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by=dx)
    knots <- round(knots, 8)
  }else{
    knots <- knots
    dx <- diff(knots)[1]
  }
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  ##
  knots1 <- knots[-c(1,length(knots))]
  P <- outer(x, knots1, MortSmooth_tpower, deg-1)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg)/(gamma(deg)*dx^(deg-1))
  BB <- ((-1)^(deg) * P %*% t(D))/dx
  D <- diff(diag(ncol(BB) + 1))
  C <- BB %*% D
  ##
  out <- list(dx=dx, knots=knots, B=B, C=C)
}


# dt_in
# # number of knots
# nd = 18
# # age interval
# ag1 = 0
# ag2 = 90
# # year
# yr = 2018
# # type of flu
# tp = "h1"

get_best_aic_bic <- function(dt_in = dt_in, 
                        # knots
                        nd = 18, 
                        # ages to include
                        ag1 = 0, 
                        ag2 = 90,
                        # year
                        yr = 2009, 
                        # measure
                        tp = "h1"){
  
  dati0 <- 
    dt_in %>% 
    rename(exposure = pop) %>% 
    filter(age %in% ag1:ag2,
           year == yr,
           type == tp)
  
  x <- dati0$age
  m <- length(x)
  y <- dati0$dts ## !!!
  e <- dati0$exposure
  lmx <- log(y/e)
  lmx1 <- emp.der(x, lmx)
  
  ## B-splines and C-matrix over age
  BC <- BsplineGrad(x, min(x), max(x), nd, 3)
  B <- BC$B
  C <- BC$C
  nb <- ncol(B)
  
  ## penalty stuff
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  
  # ## 95% CI for eta and eta1
  alpha <- qnorm(0.975)
  # 
  # ## finer grid, mainly for plotting
  ms <- 91
  xs <- seq(min(x), max(x), length=ms)
  BCs <- BsplineGrad(xs, min(x), max(x), nd, 3)
  Bs <- BCs$B
  Cs <- BCs$C
  
  ## estimating for different lambdas
  lambdas <- 10^seq(-4, 6, 0.1)
  # lambdas <- 10^seq(-4,7)
  nl <- length(lambdas)
  ## what need to be saved
  BETAS <- matrix(0, nb, nl)
  V.BETAS <- array(0, dim=c(nb,nb,nl))
  AICs <- numeric(nl)
  BICs <- numeric(nl)
  
  ## penalized IWLS algorithm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  for(l in 1:nl){
    P <- lambdas[l]*tDD
    ## penalized IWLS algorithm 
    eta <- log((y+1)/(e+1))
    for(it in 1:10){
      mu <- e*exp(eta)
      z <- eta + (y-mu)/mu
      w <- c(mu)
      tBWB <- t(B) %*% (w*B)
      tBWz <- t(B) %*% (w*z)
      tBWBpP <- tBWB+P
      betas <- solve(tBWBpP, tBWz)
      old.eta <- eta
      eta <- B %*% betas
      dif.eta <- max(abs(old.eta - eta))
      if(dif.eta < 1e-6) break
    }
    betas.hat <- betas
    V.betas <- solve(tBWBpP)
    BETAS[,l] <- betas.hat
    V.BETAS[,,l] <- V.betas
    ## BIC
    H <- V.betas%*%tBWB
    h <- diag(H)
    ed <- sum(h)
    y0 <- y
    y0[y==0] <- 10^-8
    dev <- 2 * sum((y * log(y0/mu)))
    AICs[l] <- dev + 2*ed
    BICs[l] <- dev + log(m)*ed
    cat(lambdas[l], BICs[l], it, dif.eta, "\n")
  }
  
  pmin_aic <- which.min(AICs)
  (lambda.hat_aic <- lambdas[pmin_aic])
  pmin_bic <- which.min(BICs)
  (lambda.hat_bic <- lambdas[pmin_bic])
  
  aic <- AICs[pmin_aic]
  bic <- BICs[pmin_bic]
  

  # outcome data from function 
  list_out <- list(
    aic = aic,
    bic = bic
  )
  
  return(list_out)
  
}

do_gc_magic_best <- function(dt_in = dt_in, 
                        # knots
                        nd = 18, 
                        # ages to include
                        ag1 = 0, 
                        ag2 = 90,
                        # year
                        yr = 2009, 
                        # measure
                        tp = "h1"){
  
  dati0 <- 
    dt_in %>% 
    rename(exposure = pop) %>% 
    filter(age %in% ag1:ag2,
           year == yr,
           type == tp)
  
  x <- dati0$age
  m <- length(x)
  y <- dati0$dts ## !!!
  e <- dati0$exposure
  lmx <- log(y/e)
  lmx1 <- emp.der(x, lmx)
  
  ## B-splines and C-matrix over age
  BC <- BsplineGrad(x, min(x), max(x), nd, 3)
  B <- BC$B
  C <- BC$C
  nb <- ncol(B)
  
  ## penalty stuff
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  
  # ## 95% CI for eta and eta1
  alpha <- qnorm(0.975)
  # 
  # ## finer grid, mainly for plotting
  ms <- 91
  xs <- seq(min(x), max(x), length=ms)
  BCs <- BsplineGrad(xs, min(x), max(x), nd, 3)
  Bs <- BCs$B
  Cs <- BCs$C
  
  ## estimating for different lambdas
  lambdas <- 10^seq(-4, 6, 0.1)
  # lambdas <- 10^seq(-4,7)
  nl <- length(lambdas)
  ## what need to be saved
  BETAS <- matrix(0, nb, nl)
  V.BETAS <- array(0, dim=c(nb,nb,nl))
  AICs <- numeric(nl)
  BICs <- numeric(nl)
  
  ## penalized IWLS algorithm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  for(l in 1:nl){
    P <- lambdas[l]*tDD
    ## penalized IWLS algorithm 
    eta <- log((y+1)/(e+1))
    for(it in 1:10){
      mu <- e*exp(eta)
      z <- eta + (y-mu)/mu
      w <- c(mu)
      tBWB <- t(B) %*% (w*B)
      tBWz <- t(B) %*% (w*z)
      tBWBpP <- tBWB+P
      betas <- solve(tBWBpP, tBWz)
      old.eta <- eta
      eta <- B %*% betas
      dif.eta <- max(abs(old.eta - eta))
      if(dif.eta < 1e-6) break
    }
    betas.hat <- betas
    V.betas <- solve(tBWBpP)
    BETAS[,l] <- betas.hat
    V.BETAS[,,l] <- V.betas
    ## BIC
    H <- V.betas%*%tBWB
    h <- diag(H)
    ed <- sum(h)
    y0 <- y
    y0[y==0] <- 10^-8
    dev <- 2 * sum((y * log(y0/mu)))
    AICs[l] <- dev + 2*ed
    BICs[l] <- dev + log(m)*ed
    cat(lambdas[l], BICs[l], it, dif.eta, "\n")
  }
  
  pmin_aic <- which.min(AICs)
  (lambda.hat_aic <- lambdas[pmin_aic])
  pmin_bic <- which.min(BICs)
  (lambda.hat_bic <- lambdas[pmin_bic])
  
  aic <- AICs[pmin_aic]
  bic <- BICs[pmin_bic]
  
  ## compute etas and etas1 and 95% CI for each lambda
  ETAS <- ETAS1 <- matrix(0, ms, nl)
  ETAS.UP <- ETAS.LOW <- ETAS1.UP <- ETAS1.LOW <- matrix(0, ms, nl)
  
  l <- 1
  for(l in 1:nl){
    betas.hat <- BETAS[,l]
    V.betas <- V.BETAS[,,l]
    etas.hat <- Bs%*%betas.hat
    etas1.hat <- Cs%*%betas.hat
    V.etas <- Bs %*% V.betas %*% t(Bs)
    V.etas1 <- Cs %*% V.betas %*% t(Cs)
    se.etas <- sqrt(diag(V.etas))
    se.etas1 <- sqrt(diag(V.etas1))
    etas.up <- etas.hat+alpha*se.etas
    etas.low <- etas.hat-alpha*se.etas
    etas1.up <- etas1.hat+alpha*se.etas1
    etas1.low <- etas1.hat-alpha*se.etas1
    ## saving
    ETAS[,l] <- etas.hat
    ETAS1[,l] <- etas1.hat
    ETAS.LOW[,l] <- etas.low
    ETAS.UP[,l] <- etas.up
    ETAS1.LOW[,l] <- etas1.low
    ETAS1.UP[,l] <- etas1.up
  }
  
  ## at the alpha=0.05
  ## for a given alpha we have different outcome in terms of rate-of-aging
  ## 1) not significant rate-of-aging
  ## 2) positive significant rate-of-aging
  ## 3) negative significant rate-of-aging
  
  SIGN.ETAS1 <- matrix(0, ms, nl)
  l=1
  for(l in 1:nl){
    etas1.low.l <- ETAS1.LOW[,l] 
    etas1.up.l <- ETAS1.UP[,l] 
    NoSign.l <- which(etas1.low.l<0 & etas1.up.l>0)
    Neg.l <- which(etas1.low.l<0 & etas1.up.l<0)
    Pos.l <- which(etas1.low.l>0 & etas1.up.l>0)
    SIGN.ETAS1[NoSign.l,l] <- 0
    SIGN.ETAS1[Neg.l,l] <- -1
    SIGN.ETAS1[Pos.l,l] <- 1
  }
  
  best_lambda_aic <- lambdas[pmin_aic]
  best_lambda_bic <- lambdas[pmin_bic]
  
  DFsign_aic <- 
    expand.grid(list(ages=xs, lambdas=log10(lambdas))) %>% 
    mutate(cohort = yr - ages,
           sign = c(SIGN.ETAS1) %>% as.factor) %>% 
    filter(lambdas == log10(best_lambda_aic))
  
  DFsign_bic <- 
    expand.grid(list(ages=xs, lambdas=log10(lambdas))) %>% 
    mutate(cohort = yr - ages,
           sign = c(SIGN.ETAS1) %>% as.factor) %>% 
    filter(lambdas == log10(best_lambda_bic))
  
  dt_sizer <- 
    bind_rows(
      DFsign_aic %>% 
        mutate(opt = "AIC"),
      DFsign_bic %>% 
        mutate(opt = "BIC")
    ) %>% 
    as_tibble()
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # putting all together in tidy format
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # mx_obs <- 
  #   tibble(ages=x, value=lmx) %>% 
  #   mutate(# mx = 1e3*exp(eta),
  #     type1 = "Log mx")
  # 
  # d_obs <- 
  #   tibble(ages=x, value=lmx1) %>% 
  #   mutate(type1 = "Slope")
  
  # dt_obs <- 
  #   bind_rows(mx_obs, d_obs) %>% 
  #   mutate(cohort = yr - ages, 
  #          type="Actual")
  
  spls <-
    expand_grid(ages=xs, type=lambdas) %>%
    arrange(type, ages) %>%
    mutate(value = c(ETAS),
           ll = c(ETAS.LOW),
           ul = c(ETAS.UP),
           # mx = 1e3*exp(eta),
           # mx_l = 1e3*exp(eta_l),
           # mx_u = 1e3*exp(eta_u),
           type1 = "Log mx",
           inc = 1)
  
  d_slps <-
    expand_grid(ages=xs, type=lambdas) %>%
    arrange(type, ages) %>%
    mutate(value = c(ETAS1),
           ll = c(ETAS1.LOW),
           ul = c(ETAS1.UP),
           type1 = "Slope",
           inc = ifelse(value >=-1.5 & value <=1.5, 1, 0))
  
  dt_ests <- 
    bind_rows(spls, d_slps) %>% 
    mutate(opt = case_when(type == lambda.hat_aic ~ "AIC",
                           type == lambda.hat_bic ~ "BIC",
                           TRUE ~ "no"),
           cohort = yr - ages) %>% 
    filter(opt != "no")

  # outcome data from function 
  list_out <- 
    list(
      dt_ests = dt_ests, 
      dt_sizer = dt_sizer
    )
  
  return(list_out)
  
}

# estimates for all lambdas ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do_gc_magic <- function(dt_in = dt_in, 
                        # knots
                        nd = 18, 
                        # ages to include
                        ag1 = 0, 
                        ag2 = 90,
                        # year
                        yr = 2009, 
                        # measure
                        tp = "h1"){
  
  dati0 <- 
    dt_in %>% 
    rename(exposure = pop) %>% 
    filter(age %in% ag1:ag2,
           year == yr,
           type == tp)
  
  x <- dati0$age
  m <- length(x)
  y <- dati0$dts ## !!!
  e <- dati0$exposure
  lmx <- log(y/e)
  lmx1 <- emp.der(x, lmx)
  
  ## B-splines and C-matrix over age
  BC <- BsplineGrad(x, min(x), max(x), nd, 3)
  B <- BC$B
  C <- BC$C
  nb <- ncol(B)
  
  ## penalty stuff
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  
  # ## 95% CI for eta and eta1
  alpha <- qnorm(0.975)
  # 
  # ## finer grid, mainly for plotting
  ms <- 91
  xs <- seq(min(x), max(x), length=ms)
  BCs <- BsplineGrad(xs, min(x), max(x), nd, 3)
  Bs <- BCs$B
  Cs <- BCs$C
  
  ## estimating for different lambdas
  lambdas <- 10^seq(-4, 6, 0.1)
  # lambdas <- 10^seq(-4,7)
  nl <- length(lambdas)
  ## what need to be saved
  BETAS <- matrix(0, nb, nl)
  V.BETAS <- array(0, dim=c(nb,nb,nl))
  AICs <- numeric(nl)
  BICs <- numeric(nl)
  
  ## penalized IWLS algorithm
  # ~~~~~~~~~~~~~~~~~~~~~~~~~
  for(l in 1:nl){
    P <- lambdas[l]*tDD
    ## penalized IWLS algorithm 
    eta <- log((y+1)/(e+1))
    for(it in 1:10){
      mu <- e*exp(eta)
      z <- eta + (y-mu)/mu
      w <- c(mu)
      tBWB <- t(B) %*% (w*B)
      tBWz <- t(B) %*% (w*z)
      tBWBpP <- tBWB+P
      betas <- solve(tBWBpP, tBWz)
      old.eta <- eta
      eta <- B %*% betas
      dif.eta <- max(abs(old.eta - eta))
      if(dif.eta < 1e-6) break
    }
    betas.hat <- betas
    V.betas <- solve(tBWBpP)
    BETAS[,l] <- betas.hat
    V.BETAS[,,l] <- V.betas
    ## BIC
    H <- V.betas%*%tBWB
    h <- diag(H)
    ed <- sum(h)
    y0 <- y
    y0[y==0] <- 10^-8
    dev <- 2 * sum((y * log(y0/mu)))
    AICs[l] <- dev + 2*ed
    BICs[l] <- dev + log(m)*ed
    cat(lambdas[l], BICs[l], it, dif.eta, "\n")
  }
  
  pmin_aic <- which.min(AICs)
  (lambda.hat_aic <- lambdas[pmin_aic])
  pmin_bic <- which.min(BICs)
  (lambda.hat_bic <- lambdas[pmin_bic])
  
  ## compute etas and etas1 and 95% CI for each lambda
  ETAS <- ETAS1 <- matrix(0, ms, nl)
  ETAS.UP <- ETAS.LOW <- ETAS1.UP <- ETAS1.LOW <- matrix(0, ms, nl)
  
  l <- 1
  for(l in 1:nl){
    betas.hat <- BETAS[,l]
    V.betas <- V.BETAS[,,l]
    etas.hat <- Bs%*%betas.hat
    etas1.hat <- Cs%*%betas.hat
    V.etas <- Bs %*% V.betas %*% t(Bs)
    V.etas1 <- Cs %*% V.betas %*% t(Cs)
    se.etas <- sqrt(diag(V.etas))
    se.etas1 <- sqrt(diag(V.etas1))
    etas.up <- etas.hat+alpha*se.etas
    etas.low <- etas.hat-alpha*se.etas
    etas1.up <- etas1.hat+alpha*se.etas1
    etas1.low <- etas1.hat-alpha*se.etas1
    ## saving
    ETAS[,l] <- etas.hat
    ETAS1[,l] <- etas1.hat
    ETAS.LOW[,l] <- etas.low
    ETAS.UP[,l] <- etas.up
    ETAS1.LOW[,l] <- etas1.low
    ETAS1.UP[,l] <- etas1.up
  }
  
  ## at the alpha=0.05
  ## for a given alpha we have different outcome in terms of rate-of-aging
  ## 1) not significant rate-of-aging
  ## 2) positive significant rate-of-aging
  ## 3) negative significant rate-of-aging
  
  SIGN.ETAS1 <- matrix(0, ms, nl)
  l=1
  for(l in 1:nl){
    etas1.low.l <- ETAS1.LOW[,l] 
    etas1.up.l <- ETAS1.UP[,l] 
    NoSign.l <- which(etas1.low.l<0 & etas1.up.l>0)
    Neg.l <- which(etas1.low.l<0 & etas1.up.l<0)
    Pos.l <- which(etas1.low.l>0 & etas1.up.l>0)
    SIGN.ETAS1[NoSign.l,l] <- 0
    SIGN.ETAS1[Neg.l,l] <- -1
    SIGN.ETAS1[Pos.l,l] <- 1
  }
  
  DFsign <- 
    expand.grid(list(ages=xs, lambdas=log10(lambdas))) %>% 
    mutate(cohort = yr - ages,
           sign = c(SIGN.ETAS1) %>% as.factor)
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # putting all together in tidy format
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  mx_obs <- 
    tibble(ages=x, value=lmx) %>% 
    mutate(# mx = 1e3*exp(eta),
      type1 = "Log mx")
  
  spls <- 
    expand_grid(ages=xs, type=lambdas) %>% 
    arrange(type, ages) %>% 
    mutate(value = c(ETAS),
           ll = c(ETAS.LOW),
           ul = c(ETAS.UP),
           # mx = 1e3*exp(eta),
           # mx_l = 1e3*exp(eta_l),
           # mx_u = 1e3*exp(eta_u),
           type1 = "Log mx",
           inc = 1)
  
  d_obs <- 
    tibble(ages=x, value=lmx1) %>% 
    mutate(type1 = "Slope")
  
  d_slps <- 
    expand_grid(ages=xs, type=lambdas) %>% 
    arrange(type, ages) %>% 
    mutate(value = c(ETAS1),
           ll = c(ETAS1.LOW),
           ul = c(ETAS1.UP),
           type1 = "Slope",
           inc = ifelse(value >=-1.5 & value <=1.5, 1, 0))
  
  dt_obs <- 
    bind_rows(mx_obs, d_obs) %>% 
    mutate(cohort = yr - ages, 
           type="Actual")
  
  dt_ests <- 
    bind_rows(spls, d_slps) %>% 
    mutate(opt = case_when(type == lambda.hat_aic ~ "AIC",
                           type == lambda.hat_bic ~ "BIC",
                           TRUE ~ "no"),
           cohort = yr - ages)
  
  list_out <- list(
    # p_psplines = p1,
    # p_sizer = p2,
    dt_obs = dt_obs,
    dt_ests = dt_ests,
    dt_sizer = DFsign
  )
  
  return(list_out)
  
}
