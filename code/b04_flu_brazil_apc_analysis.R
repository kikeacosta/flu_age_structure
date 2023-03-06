rm(list=ls())
{
  libs <- c("Epi", "MASS", "writexl")
  for (i in libs){
    library(i,character.only = TRUE)
  }
  library(tidyverse)
  select <- dplyr::select
}

# loading data ====
# ~~~~~~~~~~~~~~~~~.

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

sex <- "t"
sub = "h1"
num <- "cases"
den <- "population"

# size of categories
gr <- 1
# ages
amin <- 0;
amax <- 85;
# periods (maximum 2014, because the first group includes seasonal years 1959-60 & 1960-61, so
# the group includes the seasonal periods 2013-14 & 2014-15, thus, seasonal period 2015-16 is excluded)
pmin <- 2009; 
pmax <- 2019; 
# distribution
dist <- "nb"
# drift extraction (consult Crastensen material)
extr <- "n" 

parm="APC"; mod="factor";
dr.extr="n";
amn=amin; amx=amax;
pmn=pmin; pmx=pmax;
dis=dist
sb <- sub
nm <- "cases"
dn <- "population"
sx <- "t"


# function ====
## detrended models with Negative Binomial or Poisson distribution
apc_acp <- 
  function(gr = 2, 
           sb = sub,
           nm = num, 
           dn = den,
           sx = sex,
           amn = amin, 
           amx = amax, 
           pmn = pmin, 
           pmx = pmax, 
           parm = "APC", 
           dis = dist, 
           mod = "factor", 
           dr.extr = "y"){
    
    db <- 
      db_in %>% 
      filter(sex == sx,
             sub == sb,
             year >= pmn, 
             year <= pmx, 
             age >= amn, 
             age <= amx) %>% 
      mutate(D = case_when(nm == "deaths" ~ deaths,
                           nm == "cases" ~ cases,
                           nm == "population" ~ pop),
             Y = case_when(dn == "deaths" ~ deaths,
                           dn == "cases" ~ cases,
                           dn == "population" ~ pop)) %>% 
      mutate(A = as.integer(amn + floor((age - amn) / gr) * gr),
             P = as.integer(pmn + floor((year - pmn) / gr) * gr),
             D = as.integer(D),
             Y = as.integer(Y)) %>%
      select(A, P, D, Y) %>% 
      group_by(A, P) %>% 
      summarise(D = sum(D), Y = sum(Y)) %>% 
      ungroup() %>% 
      mutate(C = P - A,
             D = round(D))
    
    ############################
    # appling apc.fit function
    ############################
    
    knos <- c(round(length(unique(db$A))/5), 
              round(length(unique(db$P))/5), 
              round(length(unique(db$C))/5))
    
    alpha <- 0.05
    ref.p <- NA
    ref.c <- NA
    
    A <- db$A
    P <- db$P
    D <- db$D
    Y <- db$Y
    C <- db$P - db$A
    
    # Reference period and cohort (default values at the median value)
    med <- function(x, y) {
      o <- order(x)
      a <- y[o]
      names(a) <- x[o]
      return(as.numeric(names(a[cumsum(a)/sum(a) > 0.5][1])))
    }
    
    p0 <- ifelse(is.na(ref.p), med(P, D), ref.p)
    c0 <- ifelse(is.na(ref.c), med(P - A, D), ref.c)
    ref.p <- !is.na(ref.p)
    ref.c <- !is.na(ref.c)
    
    # function to project the columns of M on the orthogonal complement of [1|p] (a weighted inner product). Formula B1
    proj.ip <- function(X, M, orth = FALSE, weight = rep(1, nrow(X))) {
      Pp <- solve(crossprod(X * sqrt(weight)), t(X * weight)) %*% M
      PM <- X %*% Pp
      if (orth) 
        PM <- M - PM
      else PM
    }
    
    # function to shaved down to full rank (QR decomposition)
    Thin.col <- function(X, tol = 0.000001) {
      QR <- qr(X, tol = tol, LAPACK = FALSE)
      X[, QR$pivot[seq(length = QR$rank)], drop = FALSE]
    }
    
    # function for detrending the design matrix M
    detrend <- function(M, t, weight = rep(1, nrow(M))) {
      Thin.col(proj.ip(cbind(1, t), M, orth = TRUE, weight = weight))
    }
    
    # generating model matrices 
    # if (model == "factor") {
    MA <- model.matrix(~factor(A) - 1)
    MP <- model.matrix(~factor(P) - 1)
    MC <- model.matrix(~factor(C) - 1)
    Rp <- MP[abs(P - p0) == min(abs(P - p0)), , drop = FALSE][1, ]
    Rc <- MC[abs(C - c0) == min(abs(C - c0)), , drop = FALSE][1, ]
    # }
    
    #########
    # poisson regression
    #########
    
    if (dis == "poisson") {
      m.APC <- glm(D ~ MA + I(P - p0) + MP + MC, offset = log(Y),
                   family = poisson)
      # summary(m.APC)
      Dist <- "Poisson with log(Y) offset"
    }
    
    #########
    # negative binomial regression
    #########
    if (dis == "nb") {
      # m.APC.nb <- glm.nb(D ~ MA + MP + MC, offset = log(Y))
      # offset not working
      m.APC.nb <- glm.nb( D ~ factor(A) + factor(P) + factor(C) + offset( log(Y) ))
      # summary(m.APC.nb)
      m.APC <- glm(D ~ MA + I(P - p0) + MP + MC, offset = log(Y),
                   family = negative.binomial(m.APC.nb$theta))
      Dist <- "NB with log(Y) offset"
    }
    ##############
    
    m.AP <- update(m.APC, . ~ . - MC)
    m.AC <- update(m.APC, . ~ . - MP)
    m.Ad <- update(m.AP, . ~ . - MP)
    m.A <- update(m.Ad, . ~ . - I(P - p0))
    m.0 <- update(m.A, . ~ . - MA)
    AOV <- anova(m.A, m.Ad, m.AC, m.APC, m.AP, m.Ad, test = "Chisq")
    attr(AOV, "heading") <- "\\nAnalysis of deviance for Age-Period-Cohort model\\n"
    attr(AOV, "row.names") <- c("Age", "Age-drift", "Age-Cohort", 
                                "Age-Period-Cohort", "Age-Period", "Age-drift")
    
    A.pt <- unique(A)
    A.pos <- match(A.pt, A)
    P.pt <- unique(P)
    P.pos <- match(P.pt, P)
    C.pt <- unique(P - A)
    C.pos <- match(C.pt, P - A)
    
    MA <- cbind(1, MA)
    
    if (is.character(dr.extr)) {
      wt <- rep(1, length(D))
      drtyp <- "1-weights"
      if (toupper(substr(dr.extr, 1, 1)) %in% c("W", "T", "D")) {
        wt <- D
        drtyp <- "D-weights"
      }
      if (toupper(substr(dr.extr, 1, 1)) %in% c("L", "R")) {
        wt <- (Y^2)/D
        drtyp <- "Y2/D-weights"
      }
      if (toupper(substr(dr.extr, 1, 1)) %in% c("Y")) {
        wt <- Y
        drtyp <- "Y-weights"
      }
    }
    
    ##############
    # Fixing the references for period and cohort
    ##############
    Rp <- matrix(Rp, nrow = 1)
    Rc <- matrix(Rc, nrow = 1)
    xP <- detrend(rbind(Rp, MP), c(p0, P), weight = c(0, wt))
    xC <- detrend(rbind(Rc, MC), c(c0, P - A), weight = c(0, wt))
    MPr <- xP[-1, , drop = FALSE] - ref.p * xP[rep(1, nrow(MP)), 
                                               , drop = FALSE]
    MCr <- xC[-1, , drop = FALSE] - ref.c * xC[rep(1, nrow(MC)), 
                                               , drop = FALSE]
    
    # ACP / APC
    m.APC <- update(m.0, . ~ . - 1 + MA + I(P - p0) +
                      MPr + MCr)
    
    drift <- rbind(ci.exp(m.APC, subset = "I[(]", alpha = alpha), 
                   ci.exp(m.Ad, subset = "I[(]", alpha = alpha))
    
    rownames(drift) <- c(paste("APC (", drtyp, ")", sep = ""), 
                         "A-d")
    
    if (parm == "APC") {
      MPr <- cbind(P - p0, MPr)
      m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    if (parm == "ACP") {
      MCr <- cbind(P - A - c0, MCr)
      m.APC <- update(m.0, . ~ . - 1 + MA + MPr + MCr)
    }
    
    summary(m.APC)
    
    Age <- cbind(Age = A.pt, ci.exp(m.APC, subset = "MA", ctr.mat = MA[A.pos, , drop = FALSE], alpha = alpha))[order(A.pt),]
    Per <- cbind(Per = P.pt, ci.exp(m.APC, subset = "MPr", ctr.mat = MPr[P.pos, , drop = FALSE], alpha = alpha))[order(P.pt),]
    Coh <- cbind(Coh = C.pt, ci.exp(m.APC, subset = "MCr", ctr.mat = MCr[C.pos, , drop = FALSE], alpha = alpha))[order(C.pt),]
    lu <- paste(formatC(c(alpha/2, 1 - alpha/2) * 100, format = "f", 
                        digits = 1), "%", sep = "")
    
    colnames(Age)[-1] <- c("Rate", lu)
    colnames(Per)[-1] <- c("P-RR", lu)
    colnames(Coh)[-1] <- c("C-RR", lu)
    Type <- paste("ML of APC-model", Dist, ": (", parm, "):\\n")
    Model <- m.APC
    
    res <- list(Type = Type, Model = Model, Age = Age, Per = Per, 
                Coh = Coh, Drift = drift, 
                Ref = c(Per = if (parm %in% c("APC", "ADPC", "Ad-P-C", "AP-C")) p0 else NA, 
                        Coh = if (parm %in% c("ACP", "ADCP", "Ad-C-P", "AC-P")) c0 else NA), 
                Anova = AOV)
    
    age_effects <- 
      Age %>% 
      as_tibble() %>% 
      rename(value = 'Rate',
             ll = '2.5%',
             ul = '97.5%',
             Years = Age) %>% 
      group_by() %>% 
      mutate(value = value / mean(value),
             ll = ll / mean(ll),
             ul = ul / mean(ul)) %>% 
      ungroup() %>% 
      mutate(dimension = "Age")
    
    per_effects <- 
      Per %>% 
      as_tibble() %>% 
      rename(value = `P-RR`,
             ll = '2.5%',
             ul = '97.5%',
             Years = Per) %>% 
      mutate(dimension = "Period")
    
    coh_effects <- 
      Coh %>% 
      as_tibble() %>% 
      rename(value = 'C-RR',
             ll = '2.5%',
             ul = '97.5%',
             Years = Coh) %>% 
      mutate(dimension = "Cohort")
    
    apc_out <-
      bind_rows(age_effects,
                per_effects,
                coh_effects) %>% 
      mutate(model=parm,
             drift=drift[1,1],
             dimension = factor(dimension, levels = c("Age", "Period", "Cohort")))
    
    return(apc_out)
  }


# h1 effects
# ~~~~~~~~~~
{
  ## Cohort detrended model 
  apc_h1_incid <- 
    apc_acp(sx = "t", sb = "h1", 
            nm = "cases", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  apc_h1_death <- 
    apc_acp(sx = "t", sb = "h1", 
            nm = "deaths", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  apc_h3_incid <- 
    apc_acp(sx = "t", sb = "h3", 
            nm = "cases", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  apc_h3_death <- 
    apc_acp(sx = "t", sb = "h3", 
            nm = "deaths", dn = "population",
            gr = 2,
            amn = 5, 
            amx = 85, 
            parm = "APC", dr.extr = extr)
  
  
  apc_incid <- 
    apc_h1_incid %>% 
    select(Years, rr_h1 = value, dimension) %>% 
    left_join(apc_h3_incid %>% 
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

  
  ## Drift value (equal in both models)
  apc_h1_incid$drift[1] - 1
  apc_h1_death$drift[1] - 1
  apc_h3_incid$drift[1] - 1
  apc_h3_death$drift[1] - 1
  
  ## Plot of cohort effects
  apc_h1_incid %>% 
    filter(dimension == "Cohort",
           model == "APC",
           Years %in% 1930:2016) %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    geom_ribbon(aes(Years, ymin = ll, ymax = ul), 
                alpha = 0.3)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
    annotate("text", 
             y = 1, 
             x = c(1947, 1957, 1968, 1978, 2009), 
             label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
             angle = 90, hjust = 1, vjust = 0, size = 3)+
    scale_x_continuous(breaks = seq(1910, 2010, 10))+
    scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5))+
    labs(x = "Cohort", title = "H1N1")+
    theme_bw()
  
  ## Plot of cohort effects
  apc_h1_death %>% 
    filter(dimension == "Cohort",
           model == "APC",
           Years %in% 1930:2016) %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    geom_ribbon(aes(Years, ymin = ll, ymax = ul), 
                alpha = 0.3)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
    annotate("text", 
             y = 1, 
             x = c(1947, 1957, 1968, 1978, 2009), 
             label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
             angle = 90, hjust = 1, vjust = 0, size = 3)+
    scale_x_continuous(breaks = seq(1910, 2010, 10))+
    scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5))+
    labs(x = "Cohort", title = "H1N1")+
    theme_bw()
  
  coh_per_apc %>% 
    filter(model == "APC") %>% 
    ggplot()+
    geom_line(aes(Years, value))+
    geom_ribbon(aes(Years, ymin = ll, ymax = ul), 
                alpha = 0.3)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    scale_x_continuous()+
    scale_y_log10()+
    facet_grid(~dimension, scales = "free", space = "free")+
    labs(x = "Cohort", title = "H1N1")+
    theme_bw()
  
  ggsave("figures/brazil_cohort_effects_H1N1_deaths.png",
         w = 8, h = 4)
  
}


# h3 effects
# ~~~~~~~~~~
{
  ## Cohort detrended model 
  coh_per_apc <- 
    apc_acp(sb = "h3", ot = "death",
            amn = amin, amx = amax, parm = "APC", dr.extr = extr)
  
  ## Period detrended model 
  coh_per_acp <- 
    apc_acp(sb = "h3", ot = "death",
            amn = amin, amx = amax, parm = "ACP", dr.extr = extr)
  
  ## Drift value (equal in both models)
  coh_per_apc$drift[1]
  coh_per_acp$drift[1]
  
  ## Plot of cohort effects
  rr_dt <- 
    coh_per_apc %>% 
    bind_rows(coh_per_acp) %>% 
    mutate(rr = exp(Coefficient),
           rr_l = exp(ll),
           rr_u = exp(ul))
  
  rr_dt %>% 
    filter(Effect=="Cohort",
           model == "APC") %>% 
    ggplot()+
    geom_ribbon(aes(Years, ymin = rr_l, ymax = rr_u), 
                alpha = 0.3)+
    geom_line(aes(Years, rr))+
    geom_hline(yintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
    annotate("text", 
             y = 1, 
             x = c(1947, 1957, 1968, 1978, 2009), 
             label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
             angle = 90, hjust = 1, vjust = 0, size = 3)+
    scale_y_log10(breaks = c(0.3, 0.5, 1, 2, 3))+
    scale_x_continuous(breaks = seq(1910, 2010, 10))+
    labs(x = "Cohort", title = "H3N2")+
    theme_bw()
  ggsave("figures/brazil_cohort_effects_H3N2_deaths.png",
         w = 8, h = 4)
}

# b effects
# ~~~~~~~~~~
{
  ## Cohort detrended model 
  coh_per_apc <- 
    apc_acp(sb = "b", ot = "death",
            amn = amin, amx = amax, parm = "APC", dr.extr = extr)
  
  ## Period detrended model 
  coh_per_acp <- 
    apc_acp(sb = "b", ot = "death",
            amn = amin, amx = amax, parm = "ACP", dr.extr = extr)
  
  ## Drift value (equal in both models)
  coh_per_apc$drift[1]
  coh_per_acp$drift[1]
  
  ## Plot of cohort effects
  rr_dt <- 
    coh_per_apc %>% 
    bind_rows(coh_per_acp) %>% 
    mutate(rr = exp(Coefficient),
           rr_l = exp(ll),
           rr_u = exp(ul))
  
  rr_dt %>% 
    filter(Effect=="Cohort",
           model == "APC") %>% 
    ggplot()+
    geom_line(aes(Years, rr))+
    geom_ribbon(aes(Years, ymin = rr_l, ymax = rr_u), 
                alpha = 0.3)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
    annotate("text", 
             y = 1, 
             x = c(1947, 1957, 1968, 1978, 2009), 
             label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
             angle = 90, hjust = 1, vjust = 0, size = 3)+
    scale_x_continuous(breaks = seq(1910, 2010, 10))+
    scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5))+
    labs(x = "Cohort", title = "H1N1")+
    theme_bw()
  
  ggsave("figures/brazil_cohort_effects_B_deaths.png",
         w = 8, h = 4)
  
}


coeffs_dt %>% 
  filter(Effect=="Cohort") %>% 
  mutate(Coefficient = Coefficient*-1) %>% 
  ggplot()+
  geom_line(aes(Years, Coefficient, linetype=model))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = c(1947, 1957, 1968, 1978, 2009), linetype = "dashed", col = "blue")+
  annotate("text", 
           y = 0, 
           x = c(1947, 1957, 1968, 1978, 2009), 
           label = c("1947 (H1 drift)", "1957 (H2)", "1968 (H3)", "1978 (H1)", "2009 (H1)"),
           angle = 90, hjust = 1, vjust = 0, size = 3)+
  scale_x_continuous(breaks = seq(1910, 2010, 10))+
  labs(x = "Cohort", title = "H3N2")+
  theme_bw()


detach(package:MASS)
detach(package:Epi)
