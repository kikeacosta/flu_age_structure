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

rm(list = ls())
## plotting in a difference device
# options(device="X11")
## R-studio to get the same dir as for the .R file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("code/p_splines_sizer/00_functions.R")

ColorBlind  <-  c("#88CCEE", "#CC6677", "#DDCC77", "#117733", 
                  "#332288", "#AA4499", "#44AA99", "#999933", 
                  "#882255", "#661100", "#6699CC", "#888888")

ColorBlind[3]
ColorBlind[1]
ColorBlind[2]

do_magic <- function(dt_in = dt_in, yr, tp){
  
  dati0 <- 
    dt_in %>% 
    rename(exposure = pop) %>% 
    filter(year == yr,
           type == tp)
  
  x <- dati0$age
  m <- length(x)
  y <- dati0$dts ## !!!
  e <- dati0$exposure
  lmx <- log(y/e)
  lmx1 <- emp.der(x, lmx)

  ## B-splines and C-matrix over age
  nd <- floor(m/5)
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
    BICs[l] <- dev + log(m)*ed
    cat(lambdas[l], BICs[l], it, dif.eta, "\n")
  }
  
  # plot(log10(lambdas), BICs)
  pmin <- which.min(BICs)
  (lambda.hat <- lambdas[pmin])
  
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
  
  ## plotting log-mortality
  ## a subset is needed when "many" lambdas are considered
  ## construction of labels for the legend to be improved
  
  DFetas.obs <- data.frame(ages=x, type="Actual", eta=lmx)
  DFetas.hat <- expand.grid(list(ages=xs, type=lambdas))
  DFetas.hat$eta <- c(ETAS)
  DFetas <- rbind(DFetas.obs, DFetas.hat)
  DFetas$type <- factor(DFetas$type, levels = c("Actual", lambdas))
  
  loglab <- c(0.0000001, 0.0000003, 0.000001, 0.000003, 0.00001, 0.00003, 0.0001)

  DFetas$opt <- "no"
  DFetas$opt[DFetas$type==lambda.hat] <- "yes"
  DFetas$opt <- factor(DFetas$opt)
  
  ## plotting rate-of-aging
  DFetas1.obs <- data.frame(ages=x, type="Actual", eta1=lmx1)
  DFetas1.hat <- expand.grid(list(ages=xs, type=lambdas))
  DFetas1.hat$eta1 <- c(ETAS1)
  DFetas1 <- rbind(DFetas1.obs, DFetas1.hat)
  DFetas1$type <- factor(DFetas1$type, levels = c("Actual", lambdas))
  
  DFetas1$opt <- "no"
  DFetas1$opt[DFetas1$type==lambda.hat] <- "yes"
  DFetas1$opt <- factor(DFetas1$opt)
  
  ## plotting log-mortality and rate-of-aging
  ## for each lambda in the same framework
  ## Note: different y-breaks to be done
  DF1etas  <- DFetas
  DF1etas1 <- DFetas1
  names(DF1etas)[3] <- "value"
  names(DF1etas1)[3]  <- "value"
  DF1etas$type1 <- "Log-mortality"
  DF1etas1$type1 <- "Rate-of-aging"
  DFall <- 
    rbind(DF1etas, DF1etas1) %>% 
    mutate(cohort = yr - ages) 
  

  DFetas1$opt <- "no"
  DFetas1$opt[DFetas1$type==lambda.hat] <- "yes"
  DFetas1$opt <- factor(DFetas1$opt)
  
  h_line <- data.frame(
    yintercept = 0,
    type1 = "Rate-of-aging"
  )

  chts <- 
    DFall %>%
    filter(cohort == min(cohort) | cohort == max(cohort)) %>% 
    pull(cohort)
  
  coh_bks <- seq((min(chts) - min(chts)%%10), max(chts) + (10 - max(chts)%%10), 10)
    
  optims <- 
    DFall %>%
    filter(opt == "yes")
  
  obs <- 
    DFall %>%
    filter(type == "Actual")
  
  p_spl_drv <- 
    DFall %>%
    filter(type != "Actual") %>% 
    ggplot(aes(x=cohort, y=value)) +
    facet_wrap(~type1, nrow = 2, scales = "free_y")+
    geom_line(aes(group = type), alpha = 0.1)+
    geom_point(data = obs, col = "black", size = 0.8)+
    geom_line(data = optims, size = 1, col = "#ff006e", alpha = 1)+
    scale_x_continuous(breaks=coh_bks)+
    geom_hline(data = h_line, aes(yintercept=yintercept), 
               col="grey30", lty=1, lwd=0.5, linetype = "dashed")+
    geom_vline(xintercept = c(1957, 1968, 1978, 2009), linetype = "dashed",
               col = "grey40")+
    coord_cartesian(xlim = c(1920, 2020))+
    labs(x = "Cohort", y = "",
         title = paste0(tp, "_", yr))+
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.text.x  = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y  = element_text(size = 14),
          # plot.title = element_text(size=20, hjust = 0),
          strip.text = element_text(size = 14),
          strip.background = element_blank(),
          legend.position = "none",
          legend.direction = "horizontal",
          legend.key.size = unit(0.1, 'cm'), #change legend key size
          legend.key.height = unit(.2, 'cm'), #change legend key height
          legend.key.width = unit(.2, 'cm'), #change legend key width
          legend.text = element_text(size=6))
  p_spl_drv
  ggsave(paste0("figures/p_splines_sizer/ests/spsplines_deriv", 
                "_",
                tp,
                "_",
                yr,
                ".png"),
         w = 10, h = 8)
  
  DFsign <- 
    expand.grid(list(ages=xs, lambdas=log10(lambdas))) %>% 
    mutate(cohort = yr - ages,
           sign = c(SIGN.ETAS1) %>% as.factor)
    
  # DFsign$sign <- c(SIGN.ETAS1)
  # DFsign$sign <- factor(DFsign$sign)
  
  # lambdalab <- lambdas
  lambdalab <- 10^seq(-4,7)
  
  p_sz <- 
    DFsign %>%
    ggplot(aes(cohort, lambdas)) + 
    geom_tile(aes(fill=sign))+
    scale_fill_manual(name=NULL,
                      values=c("0"="#DDCC77","-1"="#88CCEE","1"="#CC6677"),
                      labels = c("Significant Negative", "Not Significant", "Significant Positive"))+
    scale_y_continuous(name="smoothing parameter (log10-scale)", expand=c(0,0),
                       breaks = log10(lambdalab), labels=lambdalab)+
    geom_hline(yintercept = log10(lambda.hat), col="black", lty=1, lwd=2)+
    geom_vline(xintercept = c(1957, 1968, 1978, 2009), linetype = "dashed",
               col = "grey40")+
    ggtitle(paste0("Rate-of-aging (SiZer map by P-splines) ", tp, "_", yr))+
    coord_cartesian(xlim = c(1920, 2020))+
    theme_bw()+
    theme(axis.title.x = element_text(size = 18),
          axis.text.x  = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          axis.text.y  = element_text(size = 14),
          legend.text = element_text(size=12),
          legend.position="top",
          plot.title = element_text(size=20, hjust = 0))
  p_sz
  ggsave(paste0("figures/p_splines_sizer/ests/psplines_sizer", 
                "_",
                tp,
                "_",
                yr,
                ".png"),
         w = 10, h = 6)
  
  list_out <- list(
    p_psplines = p_spl_drv,
    psizzer = p_sz,
    dt_ests = DFall,
    dt_sizer = DFsign
  )
  
  return(list_out)

}

dt_in <- 
  readRDS("data_inter/sample_dts_hsp_bra.rds") %>% 
  filter(age < 90)

unique(dt_in$year)

t09 <- do_magic(dt_in, 2009, "h1")
t09 <- do_magic(dt_in, 2012, "h1")
t13 <- do_magic(dt_in, 2013, "h1")
t16 <- do_magic(dt_in, 2016, "h1")
# t17 <- do_magic(dt_in, 2017, "h1")
t18 <- do_magic(dt_in, 2018, "h1")


## END
