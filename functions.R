
# Here are some functions to get started (not to literally include in the paper, just to have everything in one place)
# remotes::install_github("mpascariu/ungroup")
library(tidyverse)
library(DemoDecomp)
library(ungroup)

# ----- Helper functions ----- #
# calculate q(x) from m(x). Assumes single age intervals
mx_to_qx <- function(mx) {
  1 - exp(-mx)
}

# for sake of complemtary completeness
qx_to_mx <- function(qx){
  -log(1 - qx)
}

# calculate l(x) from q(x)
qx_to_lx <- function(qx) {
  
  # only for all-state lx
  c(1, cumprod(1 - qx))
  
}

# calculate l(x) from m(x) 
mx_to_lx <- function(mx) {
  
  # assume single ages
  c(1, exp(-cumsum(mx)))
  
}

lx_to_Lx <- function(lx){
  (lx[-1] + lx[-length(lx)]) / 2
}

#Here's an old-fashioned Sullivan function, calculates the set of for now HLE, ULE and their sum. Currently framed in terms of $m(a)$.

sully_normal <- function(mx_all, pux, type = 'h') {
  
  lx_all <- mx_to_lx(mx_all)
  Lx_all <- lx_to_Lx(lx_all)
  
  if (type == "h") {
    
    return(sum(Lx_all * (1 - pux)))  
    
  }
  
  if (type == "u") {
    
    return(sum(Lx_all * pux))
    
  }
  if (type == "t") {
    
    return(sum(Lx_all))
    
  }
  
}

# This is the rates-based version; 

sully_rates <- function(qhx, qux, phux, p0, type = "h"){
  n   <- length(qhx)
  lux <- rep(0, n + 1)
  lhx <- rep(0, n + 1)
  
  lux[1] <- p0
  lhx[1] <- 1 - p0
  for (i in 1:n){
    lux[i+1] <- lux[i] * (1 - qux[i]) + lhx[i] * phux[i]
    lhx[i+1] <- lhx[i] * (1 - qhx[i] - phux[i])
  }
  lx_all <- lhx + lux
  # simple linear avg
  Lux <- lx_to_Lx(lux)
  Lhx <- lx_to_Lx(lhx)
  pux = lux[-(n+1)] / lx_all[-(n+1)]
  Lx_all <- Lux + Lhx
  # 
  # pux <- Lux / Lx_all
  Lx_all <-  (lx_all[-1] + lx_all[1:n]) / 2
  #pux <- lux[1:n] / lx_all[1:n]
  # 
 if (type == "h") {
   
   return( sum(Lx_all * (1 - pux)))
   
 }
 
  if (type == "u") {
  
     return(sum(Lx_all * (pux)))
    
 }
 if (type == "t") {
   
   return(sum(Lx_all))
   
 }

  #sum(Lhx)
}

sully_rates2 <- function(mhx, mux, phux, p0, type = "h"){
  n   <- length(mhx)
  lux <- rep(0, n + 1)
  lhx <- rep(0, n + 1)
  
  qhx <- mx_to_qx(mhx)
  qux <- mx_to_qx(mux)
  lux[1] <- p0
  lhx[1] <- 1 - p0
  for (i in 1:n){
    lux[i+1] <- lux[i] * (1 - qux[i]) + lhx[i] * phux[i]
    lhx[i+1] <- lhx[i] * (1 - qhx[i] - phux[i])
  }
  lx_all <- lhx + lux
  # simple linear avg
  Lux <- lx_to_Lx(lux)
  Lhx <- lx_to_Lx(lhx)
  pux = lux[-(n+1)] / lx_all[-(n+1)]
  Lx_all <- Lux + Lhx
  # 
  # pux <- Lux / Lx_all
  Lx_all <-  (lx_all[-1] + lx_all[1:n]) / 2
  #pux <- lux[1:n] / lx_all[1:n]
  # 
  if (type == "h") {
    
    return( sum(Lx_all * (1 - pux)))
    
  }
  
  if (type == "u") {
    
    return(sum(Lx_all * (pux)))
    
  }
  if (type == "t") {
    
    return(sum(Lx_all))
    
  }
  
  #sum(Lhx)
}



# This gets incidence from Sullivan inputs framed in terms of $m(a)$

sully_derive_rates <- function(mx_all, pux, R_guess) {
  
  n      <- length(mx_all)
  mhx    <- mx_all / ((1 - pux) + pux * R_guess)
  mux    <- mhx * R_guess
  
  qhx    <- mx_to_qx(mhx)
  qux    <- mx_to_qx(mux)
  
  lx_all <- mx_to_lx(mx_all)
  
  # need one more element to end of pux
  pux    <- c(pux,pux[n])
  lux    <- lx_all * pux
  lhx    <- lx_all * (1 - pux)
  
  # simple closeout, as we need one more element
  # could also closeout lux, lhx with extra 0 element
  
  u_netx <- lux[2:(n+1)] - lux[1:n]
  dux    <- lux[1:n] * qux
  tux    <- u_netx + dux
  phux   <- tux / lhx[1:n]
  age    <- 1:n - 1
  
  tibble(age, qux, qhx, phux)
  
}
# Take things the other direction
rates_derive_sully <- function(qux, qhx, phux, age, p0=0){
 
  n   <- length(qhx)
  lux <- rep(0, n + 1)
  lhx <- rep(0, n + 1)
  
  lux[1] <- p0
  lhx[1] <- 1 - p0
  for (i in 1:n){
    lux[i+1] <- lux[i] * (1 - qux[i]) + lhx[i] * phux[i]
    lhx[i+1] <- lhx[i] * (1 - qhx[i] - phux[i])
  }
  lx_all <- lhx + lux
  # simple linear avg
  # Lux <- lx_to_Lx(lux)
  # Lhx <- lx_to_Lx(lhx)
  pux <- lux[-(n+1)] / lx_all[-(n+1)]
  
  mux <- qx_to_mx(qux)
  mhx <- qx_to_mx(qhx)
  mx_all <- mux * pux + mhx * (1 - pux)
  tibble(age,
         mx_all,
         pux
         )
}


# Wrap both functions to accept a single vector of arguments:

sully_rates_vec <- function(pars, type = 'h'){
  
  # dimension management, getting initial conditions
  n         <- length(pars)
  p0        <- pars[n]
  pars      <- pars[-n]
  n         <- n - 1
  dim(pars) <- c(n / 3, 3)
  # main vectors
  qhx       <- pars[, 1]
  qux       <- pars[, 2]
  phux      <- pars[, 3]
  
  sully_rates(qhx  = qhx,
              qux  = qux,
              phux = phux,
              p0   = p0,
              type = type)
}

sully_normal_vec <- function(pars, type = 'h') {
  
  n         <- length(pars)
  dim(pars) <- c(n/2,2)
  mx        <- pars[,1]
  pux       <- pars[,2]
  sully_normal(mx_all = mx,
               pux    = pux,
               type   = type)
}


# A full service rates decomp function

sully_rates_decomp <- function(mx1, mx2, pux1, pux2, R1, R2, type = 'h') {
  
  init1 <- pux1[1]
  init2 <- pux2[1]
  p1    <- sully_derive_rates(mx1,pux1,R1)
  p2    <- sully_derive_rates(mx2,pux2,R2)
  cc    <- horiuchi(sully_rates_vec,
                    c(p1$qhx,
                      p1$qux,
                      p1$phux,
                      init1),
                    c(p2$qhx,
                      p2$qux,
                      p2$phux,
                      init2),
                    type = type,
                    N=20)  
  n       <- length(cc)
  init    <- cc[n]
  cc      <- cc[-n]
  n       <- n - 1
  dim(cc) <- c(n / 3, 3)  
  colnames(cc) <- c("qhx","qux","onset")
  age <- 1:nrow(cc)-1
  out <- as_tibble(cc) |> 
    mutate(age = age, .before = 1) |> 
    pivot_longer(-age, 
                 names_to  = "component", 
                 values_to = "value") |> 
    add_row(age       = 0,
            component = "init", 
            value     = init)
  
  return(out)
}


# A full service normal decomp function for an old-fashioned Sullivan function, calculates the set of for now HLE, ULE and their sum.

sully_normal_decomp <- function(mx1, mx2, pux1, pux2, type = 'h') {
  
  n  <- length(mx1)
  cc <- horiuchi(sully_normal_vec,
                 c(mx1, pux1),
                 c(mx2, pux2),
                 type = type,
                 N    = 20)
  dim(cc)      <- c(n, 2)
  colnames(cc) <- c("mx", "pux")
  age          <- 1:nrow(cc) - 1
  out          <- as_tibble(cc) |>
    mutate(age = age) |>
    pivot_longer(-age, 
                 names_to  = "component", 
                 values_to = "value")
  
  return(out)
  
}

