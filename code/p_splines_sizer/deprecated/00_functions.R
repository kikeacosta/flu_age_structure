rm(list=ls())
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
