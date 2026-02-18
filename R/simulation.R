# ---------------------------------------------------------------------------- #
library(tidyverse)
library(DemoDecomp)
library(expm)
# ---------------------------------------------------------------------------- #
# Conceptual assumption:
# Differences between populations (or time points) are assumed to be driven
# primarily by shifts in the level of the hazards (intercepts), rather than
# by changes in the age-gradient (slopes).
# We treat most improvements/deteriorations as vertical shifts
# in the hazard schedule rather than structural age-pattern changes.

# Incidence go up with age
# Healthy mortality go up and is Gompertz-like
# Recovery goes down with age
# Unhealthy mortality > healthy mortality

# Interpretation of intercepts (baseline hazard at age 50 midpoint):
# exp(intercept) gives the approximate hazard at age 50.
# These values were chosen to produce plausible age-specific dynamics:
# - Incidence increases with age
# - Healthy mortality is lower than unhealthy mortality
# - Recovery declines with age (negative slope)

# Intercepts
# exp(-3.5) = 0.030   -> 3% incidence
# exp(-7.5) = 0.00055 -> low mortality healthy
# exp(-1.5) = 0.22    -> strong recovery at 50
# exp(-6.5) = 0.0015  -> higher mortality if unhealthy

# Negative slope
# slope     = -0.06   -> recovery goes down with age

# Abbreviations
# Log hazards scale
# h - healthy
# u - unhealthy
# d - dead

# ---------------------------------------------------------------------------- #
# Initial parameters for log hazards
ground_pars <- tibble(
  trans     = c("hu", "hd", "uh", "ud"),
  intercept = c(-3.5, -7.5, -1.5, -4.5),
  slope     = c(0.07, 0.10, -0.06, 0.08)
)

# age range
age <- 50:100
# ---------------------------------------------------------------------------- #
# Expand hazards from intercepts and slopes
expand_hazards <- function(pars = ground_pars, 
                           age = 50:100,
                           age_int = 1) {
  
  expand_grid(
    age   = age,
    trans = pars$trans) %>%
    left_join(pars, by = "trans") %>%
    mutate(
      # Gompertz
      # predict for the midpoints of the interval       
      log_hazard = intercept + slope * (age - min(age) + age_int / 2), # + half of the interval
      hazard     = exp(log_hazard)
    ) %>%
    select(age, trans, hazard)
}

# Here we have hazards
# age interval is for midpoint
ground_haz <- expand_hazards(ground_pars,
                             age = 50:100,
                             age_int = 1)

# Log hazards
# Visual check. Looks good by me.
ground_haz %>%
  mutate(trans = str_to_upper(trans),
         From  = str_sub(trans, 1, 1),
         To    = str_sub(trans, 2, 2)) %>%
  ggplot(aes(x = age, y = hazard, color = To)) +
  geom_line() + 
  scale_y_log10() +
  facet_wrap( ~ From) + 
  theme_bw() + 
  theme(legend.position = "bottom")

# Mortality risk ratio
Ra <- ground_haz %>%
  filter(trans %in% c("ud", "hd")) %>% 
  pivot_wider(names_from  = trans,
              values_from = hazard) %>% 
  mutate(Ra = ud / hd) %>% 
  select(age, Ra)

# risk ratio declines with age.
Ra %>% 
  ggplot(aes(x = age, y = Ra)) + 
  geom_line() + 
  theme_bw()
# ---------------------------------------------------------------------------- #
# Rates to probabilities
# Q matrix construction
make_Q <- function(x) {
  with(as.list(x), 
       matrix(c(
         # diagonals are - sum row
         -(hu + hd),  hu,  hd,
         uh, -(uh + ud), ud,
         0,   0,   0
         ), nrow = 3, byrow = TRUE)
  )
}
# Convert hazards to transition probability matrices
haz_to_probs <- function(hazard, age) {
  
  out <- hazard %>% 
    group_nest(age) %>%
    mutate(Q = map(data, ~ make_Q(deframe(.x))),
           P = map(Q, ~ .x %>%
                     expm() %>%
                     as.data.frame() %>%
                     set_names(c("H", "U", "D")) %>%
                     cbind("from" = c("H", "U", "D")))) %>%
    select(age, P) %>%
    unnest(P)
  
  return(out)
}

# resulting probabilities
ground_P <- haz_to_probs(ground_haz, age)

# visual check. Looks fine by me.
ground_P %>% 
  pivot_longer(-c(age, from),
               names_to  = "to",
               values_to = "p") %>% 
  filter(from != "D") %>% 
  ggplot(aes(x = age, y = p, color = to)) +
  geom_line() + 
  facet_wrap(~ from) + 
  theme_bw() + 
  theme(legend.position = "bottom")
# ---------------------------------------------------------------------------- #
# Construct MSLT
# initial lh and lu
init_constant <- function(x) {
  x    <- unlist(x[c("HH", "UH", "HU", "UU")]) 
  u    <- matrix(x, nrow = 2, byrow = TRUE)
  v    <- eigen(u)$vectors[, 1]
  init <- v / sum(v)
  setNames(init, c("H", "U"))
}

# init <- c(H=0.9,U=0.1)
# Create multistate lifetable
calculate_lt <- function(P, init = NULL) {
  
  P %>% 
    pivot_longer(-c(age, from),
                 names_to  = "to",
                 values_to = "p") %>% 
    filter(from != "D") %>% 
    unite("state", c("from", "to"), sep = "") %>%
    pivot_wider(names_from  = state,
                values_from = p) %>%
    group_modify( ~ {
      .x$lh <- numeric(nrow(.x))
      .x$lu <- numeric(nrow(.x))
      
      # Initialize
      if (is.null(init)){
         init <- init_constant(.x[1, ])
      }
      .x$lh[1] <- init["H"]
      .x$lu[1] <- init["U"]
      
      # Recursive loop with initial pars
      for (i in seq_len(nrow(.x) - 1)) {
        .x$lh[i + 1] <- .x$lh[i] * .x$HH[i] + .x$lu[i] * .x$UH[i]
        .x$lu[i + 1] <- .x$lu[i] * .x$UU[i] + .x$lh[i] * .x$HU[i]
      }
      .x # return
    }) %>% 
    # LT quantities
    mutate(Ra         = UD / HD,
           prevalence = lu / (lh + lu),
           lx         = lh + lu,
           Lx         = (lx + lead(lx, default = 0)) / 2)
  }

# result  
lt <- calculate_lt(ground_P)

# survival looks fine by me
lt %>% 
  select(age, lh, lu, lx) %>% 
  pivot_longer(-c(age),
               names_to  = "var",
               values_to = "val") %>%
  ggplot(aes(x = age, y = val, color = var)) +
  geom_line() +
  theme_bw()

# life values expectancy look fine
lt %>%
  summarise(HLE = sum(lh),
            ULE = sum(lu),
            LE  = sum(lx))

# prevalence looks ok too
lt %>%
  select(age, prevalence) %>%
  ggplot(aes(x = age, y = prevalence)) +
  geom_line() +
  theme_bw()
# ---------------------------------------------------------------------------- #
# combine all steps to create mslt with one function
run_mslt <- function(pars_df, age, init = NULL) {
  
  haz <- expand_hazards(pars_df, age, age_int = 1)
  P   <- haz_to_probs(haz, age)
  lt  <- calculate_lt(P, init = init)
  
  return(lt)
}
# ---------------------------------------------------------------------------- #
# create ground summary lx and prevalence measures for optim
ground_summary <- run_mslt(ground_pars, age, init = NULL) %>%
  select(age, lx, prevalence)
# ---------------------------------------------------------------------------- #
# Run optimization for lx and prevalence
# optimize only 3 intercepts, keep slopes intact
mslt_min_summary <- function(par, 
                             fixed_trans = "uh",
                             fixed_value = -10,
                             ground_summary,
                             age  = 50:100,
                             init = NULL) {
  # initial pars
  pars <- tibble(
    trans     = c("hu", "hd", "uh", "ud"),
    intercept = c(-3.5, -7.5, -1.5, -4.5),
    slope     = c(0.07, 0.10, -0.06, 0.08)
  )
  
  # set fixed intercept
  pars$intercept[pars$trans == fixed_trans] <- fixed_value
  
  # optimize remaining 3
  free_trans <- setdiff(pars$trans, fixed_trans)
  pars$intercept[match(free_trans, pars$trans)] <- par
  
  # run model
  lt_mod <- run_mslt(pars_df = pars, age = age, init = init) %>% 
    select(age, lx1 = lx, prevalence1 = prevalence) %>% 
    left_join(ground_summary, by = "age") %>% 
    mutate(loss_lx   = (lx1 - lx) ^ 2,
           loss_prev = (prevalence1 - prevalence) ^ 2)
  
  # return single number to minimize
  return(sum(lt_mod$loss_lx + lt_mod$loss_prev))
  
}

# run optim with BFGS
result <- optim(
  par            = ground_pars$intercept[-3],  # starting guesses for hu, hd, ud
  fn             = mslt_min_summary,
  fixed_trans    = "uh",
  fixed_value    = -5,
  ground_summary = ground_summary,
  method         = "BFGS",
  control        = list(maxit = 2000)
)

# check convergence. Converged
result$convergence
result$message

# The loss is not zero because we attempt to reproduce 102 age-specific
# quantities (lx and prevalence at each age) using only a small number
# of parameters. The system is therefore overdetermined.
# Even after optimization, an exact match is not achieved.
result$value
# ---------------------------------------------------------------------------- #
# reconstruct new parameter tibble
reconstruct_pars <- function(result,
                             fixed_trans = "uh",
                             fixed_value = -5) {
  
  pars <- tibble(
    trans     = c("hu", "hd", "uh", "ud"),
    intercept = c(-3.5, -7.5, -1.5, -4.5),
    slope     = c(0.07, 0.10, -0.06, 0.08)
  )
  
  # set fixed intercept
  pars$intercept[pars$trans == fixed_trans] <- fixed_value
  
  # fill optimized intercepts
  free_trans <- setdiff(pars$trans, fixed_trans)
  pars$intercept[match(free_trans, pars$trans)] <- result$par
  
  return(pars)
}
# ---------------------------------------------------------------------------- #
# here are new parameters
new_pars <- reconstruct_pars(result,
                             fixed_trans = "uh",
                             fixed_value = -5)

# lets compare lifetables
ground_lt <- run_mslt(ground_pars, age, init = NULL)
new_lt    <- run_mslt(new_pars,    age, init = NULL)

# lx
plot(ground_lt$lx)
lines(new_lt$lx)

# prevalence
plot(ground_lt$prevalence)
lines(new_lt$prevalence)

# sums are similar, but not exact due to age identifiability
# too many values to be modeled with too little parameters
sum(ground_lt$lx) + sum(ground_lt$prevalence)
sum(new_lt$lx)    + sum(new_lt$prevalence)