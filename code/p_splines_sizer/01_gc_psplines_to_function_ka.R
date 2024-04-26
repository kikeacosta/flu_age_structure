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
    ag2 = 100,
    # year
    2018,
    # type of flu
    "h1")

# plot of splines and slopes
t09$p_psplines
# plot of splines and slopes
t09$p_sizer

# # estimates of p-splines (type1 indicates whether are Log mx
# estimates of slopes)
# t09$dt_ests


# # data for SiZer (ALAIN)
data_sizer<-t09$dt_sizer

sizer <-
  data_sizer %>%
  mutate(slope_d1 = case_when(sign ==  1 ~ "1-pos",
                              sign == -1 ~ "2-neg",
                              TRUE ~ "3-flat" ))

pal_color <- c("red","blue","purple")

sizer %>%
  ggplot()+
  geom_tile(aes(cohort, lambdas, fill = slope_d1))+
  scale_fill_manual(values = pal_color, limits = names(pal_color))+
  geom_hline(yintercept= c(t09[["lambdaminAIC"]]), col="white", linetype = "dotted")+
  geom_hline(yintercept= c(t09[["lambdaminBIC"]]), col="white", linetype = "dotted")+
  geom_vline(xintercept = c(1918,1957,1968,1978, 2009))+
  geom_vline(xintercept = c(1933,1985,1995,2003),
             linetype = "dashed")+
  annotate("text", x = 2010, y = t09[["lambdaminAIC"]],
           label = "AIC",col="white", vjust = -.2)+
  annotate("text", x = 2010, y = t09[["lambdaminBIC"]],
             label = "BIC", col = "white", vjust = -.2)+
  scale_x_reverse(breaks = seq(1900, 2020, 10), limits=c(2025, 1910))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title="Sign of 1st derivative")+
  ylab(expression(lambda))+
  coord_cartesian(expand = 0)



# testing different combinations of knots and lambdas

k <- seq(1,40,1)
minaic <- vector(mode = "numeric", length = length(k) )
minbic <- vector(mode = "numeric", length = length(k) )
lambdaminaic <- vector(mode = "numeric", length = length(k) )
lambdaminbic <- vector(mode = "numeric", length = length(k) )


annee <-2018
sub <- "h1"
a1 <- 0
a2 <- 90
  
for(kn in k){
  temp<-do_gc_magic(dt_in, kn, ag1 = a1, ag2 = a2, annee, sub)
  minaic[kn] <- temp[["minAIC"]]
  minbic[kn] <- temp[["minBIC"]]
  lambdaminaic[kn] <- temp[["lambdaminAIC"]] 
  lambdaminbic[kn] <- temp[["lambdaminBIC"]] 
}


kaicbic <- data.frame(k,minaic,lambdaminaic,minbic,lambdaminbic)

kaic <- 
  kaicbic %>% 
  arrange(minaic)

kbic <- 
  kaicbic %>% 
  arrange(minbic) %>% 
  mutate(minaic=round(minaic,2),
         minbic=round(minbic,2))

kbic$k[1]
kbic$minbic[1]
kbic$k[2]
kbic$minbic[2]


#plot minaic and minbic by k
plot_kaicbic <-
  gather(kaicbic,estimator,value,minaic,minbic)

plot_kaicbic %>% 
  ggplot()+
  geom_line(aes(k,value,col=estimator))+
  scale_x_continuous(breaks = seq(0, 90, 5))


#plot 4 sizer for best bic
n <- seq(1,5,1)
for(kn in n){
  bestbic<-do_gc_magic(dt_in, kbic$k[kn], ag1 = a1, ag2 = a2, annee, sub)
  data_sizer<-bestbic$dt_sizer

  sizer <-
    data_sizer %>%
    mutate(Sign = case_when(sign ==  1 ~ "positive",
                              sign == -1 ~ "negative",
                              TRUE ~ "zero" ))

  pal_color <- c("blue","red","purple")

  sizer %>%
    ggplot()+
    geom_tile(aes(cohort, lambdas, fill = Sign))+
    scale_fill_manual(values = pal_color, limits = names(pal_color))+
    geom_hline(yintercept= c(bestbic[["lambdaminAIC"]]), col="white", linetype = "dotted")+
    geom_hline(yintercept= c(bestbic[["lambdaminBIC"]]), col="white", linetype = "dotted")+
    geom_vline(xintercept = c(1918,1957,1968,1978, 2009))+
    geom_vline(xintercept = c(1933,1985,1995,2003),
             linetype = "dashed")+
    annotate("text", x = 2010, y = bestbic[["lambdaminAIC"]], 
           label = "",col="white", vjust = -.2)+
    annotate("text", x = 2000, y = bestbic[["lambdaminBIC"]], 
           label = "Best BIC", col = "white", vjust = -.2)+
    scale_x_reverse(breaks = seq(1900, 2020, 10), limits=c(2025, 1910))+ 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(title= paste0("Sign of the first derivative"),
       subtitle = paste0("Cases, ", bestbic[["tp"]], ", ", annee, ", ", bestbic[["nd"]], " knots, Ages ", 
                         a1, "-", a2, ", AIC=", kbic$minaic[kn],
                         ", BIC=", kbic$minbic[kn]))+
    labs(y = expression(paste(Log[10] , " (", lambda, ")")))+
    coord_cartesian(expand = 0)


  ggsave(paste0("figures/Best/Sizer_best_bic",
                "_dts_",                      #################
                sub,
                "_",
                annee,
                "_Model_",
                kn,
                "_knts_",
                bestbic[["nd"]],
                "_ages_",
                a1, "_", a2,
                ".png"),
         w = 10, h = 6)
}


## END




