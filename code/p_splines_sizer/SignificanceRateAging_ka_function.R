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
## plotting in a difference device
# options(device="X11")
## R-studio to get the same dir as for the .R file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("code/p_splines_sizer/00_functions.R")

dt_in <- 
  readRDS("data_inter/sample_dts_hsp_bra.rds")

nd = 18
# ages to include
ag1 = 0
ag2 = 90
# year
yr = 2016
# measure
tp = "h1"

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
  ms <- 500
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
  
  # plotting 
  # ~~~~~~~~
  
  x_bks <- 
    tibble(coh_bks = rev(seq(1920, 2020, 10))) %>% 
    mutate(age_bks = yr - coh_bks,
           x_bks = paste0(coh_bks, "\n(", age_bks, ")")) %>% 
    pull(x_bks)
  
  
  # p-splines and 1st derivatives
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  h_line <- data.frame(
    yintercept = 0,
    type1 = "Slope"
  )
  
  p1 <- 
    dt_ests %>%
    filter(inc == 1) %>% 
    ggplot(aes(x=cohort, y=value)) +
    facet_wrap(~type1, nrow = 2, scales = "free_y")+
    geom_ribbon(data = .  %>% filter(opt %in% c("AIC", "BIC")), 
                aes(ymin = ll, ymax = ul, fill = opt), alpha = 0.1)+
    geom_line(aes(group = type), alpha = 0.1)+
    geom_point(data = dt_obs, col = "black", size = 0.8)+
    geom_line(data = .  %>% filter(opt %in% c("AIC", "BIC")), 
              aes(col = opt), size = .8, alpha = .9)+
    scale_x_reverse(breaks = rev(seq(1920, 2020, 10)), labels = x_bks)+
    scale_fill_manual(values = c("#ff006e", "#8338ec"))+
    scale_color_manual(values = c("#ff006e", "#8338ec"))+
    geom_hline(data = h_line, aes(yintercept=yintercept), col="grey30", lty=1, lwd=0.5)+
    geom_vline(xintercept = c(1957, 1968, 1978, 2009), lty = "dashed",
               col = "grey40")+
    coord_cartesian(xlim = c(2020, 1920))+
    labs(x = "Cohort\n(Age)", 
         title = paste0(tp, "_", yr, "_knts_", nd,
                        ", ages ", ag1, " to ", ag2),
         color = "Best", fill = "Best")+
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text.x  = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y  = element_text(size = 12),
          plot.title = element_text(size=20, hjust = 0),
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(.2, 'cm'), #change legend key height
          legend.key.width = unit(.6, 'cm'), #change legend key width
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
  p1
  
  ggsave(paste0("figures/p_splines_sizer/ests/spsplines_deriv", 
                "_",
                tp,
                "_",
                yr,
                "_knts_",
                nd,
                "_ages_",
                ag1, "_", ag2,
                ".png"),
         w = 10, h = 8)
  
  lambdalab <- 10^seq(-4,7)
  
  p2 <- 
    DFsign %>%
    ggplot(aes(cohort, lambdas)) + 
    geom_tile(aes(fill=sign))+
    scale_fill_manual(name=NULL,
                      values=c("0"="#DDCC77","-1"="#88CCEE","1"="#CC6677"),
                      labels = c("Significant Negative", "Not Significant", "Significant Positive"))+
    scale_y_continuous(name="smoothing parameter (log10-scale)", expand=c(0,0),
                       breaks = log10(lambdalab), labels=lambdalab)+
    scale_x_reverse(breaks = rev(seq(1920, 2020, 10)), labels = x_bks)+
    geom_hline(yintercept = log10(lambda.hat_aic), col="#ff006e", lty=1, lwd=1)+
    geom_hline(yintercept = log10(lambda.hat_bic), col="#8338ec", lty=1, lwd=1)+
    geom_vline(xintercept = c(1957, 1968, 1978, 2009), linetype = "dashed",
               col = "grey40")+
    annotate("text", x = 2020, y = log10(lambda.hat_aic), 
             label = "AIC", col = "#ff006e", vjust = -.2)+
    annotate("text", x = 2020, y = log10(lambda.hat_bic), 
             label = "BIC", col = "#8338ec", vjust = -.2)+
    labs(title = paste0("Slope (SiZer map by P-splines)"),
         subtitle = paste0(tp, ", ", yr, ", with ", nd, " knots. Ages ", ag1, "-", ag2))+
    coord_cartesian(xlim = c(2020, 1920))+
    theme_bw()+
    theme(axis.title.x = element_text(size = 18),
          axis.text.x  = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y  = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position="top",
          plot.title = element_text(size=18, hjust = 0),
          plot.subtitle = element_text(size=16, hjust = 0),
    )
  p2
  ggsave(paste0("figures/p_splines_sizer/ests/psplines_sizer", 
                "_",
                tp,
                "_",
                yr,
                "_knts_",
                nd,
                "_ages_",
                ag1, "_", ag2,
                ".png"),
         w = 10, h = 6)
  
  list_out <- list(
    p_psplines = p1,
    psizzer = p2,
    dt_obs = dt_obs,
    dt_ests = dt_ests,
    dt_sizer = DFsign
  )
  
  return(list_out)

}

t16 <- do_gc_magic(dt_in, 5, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 10, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 15, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 18, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 19, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 30, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 50, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 75, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 100, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 200, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 500, ag1 = 0, ag2 = 90, 2016, "h1")

t16 <- do_gc_magic(dt_in, 18, ag1 = 0, ag2 = 90, 2016, "h1")
t16 <- do_gc_magic(dt_in, 18, ag1 = 0, ag2 = 110, 2016, "h1")
t16 <- do_gc_magic(dt_in, 18, ag1 = 20, ag2 = 90, 2016, "h1")



unique(dt_in$year)

t09 <- do_gc_magic(dt_in, 30, ag1 = 0, ag2 = 90, 2009, "h1")
t12 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2012, "h1")
t13 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2013, "h1")
t16 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2016, "h1")
t17 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2017, "h1")
t18 <- do_gc_magic(dt_in, 20, ag1 = 0, ag2 = 90, 2018, "h1")

# looking at specific outputs
t09$dt_ests
t09$dt_sizer
t09$p_psplines
t09$psizzer


## END
