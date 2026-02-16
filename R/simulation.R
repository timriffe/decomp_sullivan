# ---------------------------------------------------------------------------- #
library(tidyverse)
library(DemoDecomp)
library(expm)
# ---------------------------------------------------------------------------- #
# Assumptions: changes should come from intercept rather than slope
# 10% from slope and 90% from intercept
# change that time 2 is doing better than time 1
# change to only 1 time. And run it 2 times

# Incidence go up with age
# Healthy mortality go up and is Gompertz-like
# Recovery goes down with age
# Unhealthy mortality > healthy mortality

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

# init <- c(H=1,U=0)
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
run_mslt <- function(pars_df, age) {
  
  haz <- expand_hazards(pars_df, age, age_int = 1)
  P   <- haz_to_probs(haz, age)
  lt  <- calculate_lt(P)
  
  lt
}
# ---------------------------------------------------------------------------- #
# create ground summary LE measures for optim
ground_summary <- lt %>%
  summarise(
    HLE = sum(lh),
    ULE = sum(lu),
    LE  = sum(lx)
  )
# ---------------------------------------------------------------------------- #
# this function just takes the corresponding parameter
# as a starting parameter for optimization
# Ex: if the world is no recovery, then intercepts
# are from u = -3.5, ud =  -4.5
# I currently work with intercepts only
get_start <- function(world) {
  
  if (world == "no_recovery") {
    ground_pars %>%
      filter(trans %in% c("hu","ud")) %>%
      pull(intercept)
  } else if (world == "high_incidence") {
    ground_pars %>%
      filter(trans %in% c("hu","hd")) %>%
      pull(intercept)
  } else if (world == "mortality_ratio_shift") {
    ground_pars %>%
      filter(trans %in% c("ud","hd")) %>%
      pull(intercept)
  } else if (world == "incidence_recovery_tradeoff") {
    ground_pars %>%
      filter(trans %in% c("hu","uh")) %>%
      pull(intercept)
  }
}
# ---------------------------------------------------------------------------- #
# define our worlds names
worlds <- c("no_recovery",
            "high_incidence",
            "mortality_ratio_shift",
            "incidence_recovery_tradeoff")

# function that changes our pars with accordance to worlds
# e.g world = no recovery, then we start from the
# uh = -20, slope = 0
rebuild_pars <- function(par_vec, world) {
  
  pars_mod <- ground_pars
  
  if (world == "no_recovery") {
    
    # remove recovery properly
    pars_mod$intercept[pars_mod$trans == "uh"] <- -20
    pars_mod$slope[pars_mod$trans == "uh"]     <- 0
    
    # optimize unhealthy mortality
    pars_mod$intercept[pars_mod$trans == "ud"] <- par_vec[1]
    pars_mod$slope[pars_mod$trans == "ud"]     <- par_vec[2]
  }
  
  if (world == "high_incidence") {
    
    # optimize incidence
    pars_mod$intercept[pars_mod$trans == "hu"] <- par_vec[1]
    pars_mod$slope[pars_mod$trans == "hu"]     <- par_vec[2]
  }
  
  if (world == "mortality_ratio_shift") {
    
    # vary unhealthy mortality only
    pars_mod$intercept[pars_mod$trans == "ud"] <- par_vec[1]
    pars_mod$slope[pars_mod$trans == "ud"]     <- par_vec[2]
  }
  
  if (world == "incidence_recovery_tradeoff") {
    
    pars_mod$intercept[pars_mod$trans == "hu"] <- par_vec[1]
    pars_mod$slope[pars_mod$trans == "uh"]     <- par_vec[2]
  }
  
  pars_mod
}
# ---------------------------------------------------------------------------- #
# Optimization function
# I try to find the different trajectories
# that result in same HLE and LE values
# optimize the HLE and LE difference
mslt_min_summary <- function(par_vec,
                             world,
                             ground_summary,
                             age = 50:100) {
  
  pars_mod <- rebuild_pars(par_vec, world)
  
  lt_mod <- run_mslt(pars_mod, age)
  
  summ_mod <- lt_mod %>%
    summarise(
      HLE = sum(lh),
      LE  = sum(lx)
    )
  
  loss <- (summ_mod$HLE - ground_summary$HLE)^2 +
    (summ_mod$LE  - ground_summary$LE)^2
  
  return(loss)
}

# run optimization
results <- map(worlds, function(w) {
  
  optim(
    par     = get_start(w),
    fn      = mslt_min_summary,
    world   = w,
    ground_summary = ground_summary,
    method  = "Nelder-Mead",
    control = list(maxit = 2000)
  )
  
}) %>%
  set_names(worlds)

# so, currently the problem is the no_recovery 
# it did not converge, I think that setting recovery to 0
# is oo much, too strong of a constraing
# other worlds did good
rbind(map(results, "convergence"))
rbind(map(results, "value"))
rbind(map(results, "message"))

# obtained parameters
map_dfr(results, "par")
# ---------------------------------------------------------------------------- #
# summary stats
# exactly we can see that no_recovery did not converge 
# while others did
map_dfr(names(results), function(w){
  
  pars_mod <- rebuild_pars(results[[w]]$par, w)
  lt_mod   <- run_mslt(pars_mod, age)
  
  lt_mod %>%
    summarise(world = w,
              HLE   = sum(lh),
              ULE   = sum(lu),
              LE    = sum(lx))
})

# rebuild optimized parameter tables
optimized_pars <- map2(results,
                       names(results),
                       ~ rebuild_pars(.x$par, .y))

# run mslt for all worlds with new pars
lt_worlds <- map2_dfr(optimized_pars,
                      names(optimized_pars),
                      ~ run_mslt(.x, age) %>%
                        mutate(world = .y))

# check prevalence
ggplot(lt_worlds, aes(x = age, y = prevalence, color = world)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(title = "Prevalence by World") + 
  theme(legend.position = "bottom")

# mortality risk ratios
lt_worlds %>% 
  ggplot(aes(x = age, y = Ra, color = world)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(title = "Mortality Risk Ratio (Ra) by World") + 
  theme(legend.position = "bottom")
