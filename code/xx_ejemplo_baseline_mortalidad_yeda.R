rm(list=ls())
library(tidyverse)

# Mortality data from Brazil for 
dt <- read_rds("data_inter/brazil_monthly_master.rds")


# seleccionando muertes por todas klas causas entre 2015 y 2020, 
# para sexo total y edad 60
dt2 <- 
  dt %>% 
  filter(year %in% 2015:2020,
         sex == "t",
         age == 60,
         cause == "total") 
  
dt2 %>% 
  ggplot()+
  geom_line(aes(date, dts))

# agregando de mensual a anual
dt3 <- 
  dt2 %>% 
  group_by(year) %>% 
  summarise(dts = sum(dts),
            pop = mean(pop)) %>% 
  ungroup() %>% 
  mutate(mx = dts/pop,
         # weights (para decir al modelo que periodios incluir en el modelo)
         # en este caso, queremos predecir 2020, entonces lo excluimos del 
         # ajuste
         w = ifelse(year < 2020, 1, 0))


# plot de muertes anuales
dt3 %>% 
  ggplot()+
  geom_line(aes(year, dts))

# plot de tasas anuales
dt3 %>% 
  ggplot()+
  geom_line(aes(year, mx))


# modelo de regresion lineal para predecir tasas esperadas en 2015

mod <- 
  lm(mx ~ year,
     # le indicamos que tenemos weights
     weights = w,
     data = dt3)

res <- 
  predict(mod,
          newdata = dt3,
          se.fit = TRUE)

# valores ajustados (tasas de mortalidad esperadas en cada año)
res$fit
# desviacion estandard de cada valor predicho
res$se.fit

# asignando predicciones e intervalos de confianza a nuestros datos
dt4 <- 
  dt3 %>% 
  mutate(pred = res$fit,
         ul = res$fit + 1.96*res$se.fit,
         ll = res$fit - 1.96*res$se.fit)

# plot de la mortalidad observada y esperada
dt4 %>% 
  ggplot()+
  geom_ribbon(aes(year, ymin = ll, ymax = ul), alpha = 0.2)+
  geom_point(aes(year, mx))+
  geom_line(aes(year, pred))


# exceso de mortalidad:
# tasa de mortalidad observada en 2020
obs <- dt4 %>% filter(year == 2020) %>% pull(mx)
# tasa esperada en 2020
esp <- dt4 %>% filter(year == 2020) %>% pull(pred)
# población en 2020 
pop <- dt4 %>% filter(year == 2020) %>% pull(pop)

# exceso en tasa = observada - esperada
exc_mx <- obs-esp

# exceso en numero
exc_dts <- exc_mx*pop

exc_dts
