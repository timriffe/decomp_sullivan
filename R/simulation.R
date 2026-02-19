# ------------------------------------------------------------------------------
#
# Goal:
#   Build a small set (5-8) of multistate "worlds" that match the same Sullivan
#   outputs (lx, prevalence), including ~2-3 "no-return" worlds (UH ~ 0).
#
# Method:
#   (1) Define a ground world.
#   (2) Compute its Sullivan outputs.
#   (3) Find near-null directions (SVD of Jacobian of (lx, prevalence)).
#   (4) Create candidate worlds by stepping in (dir1, dir2).
#   (5) Polish each candidate using ridge-regularized bounded optimization.
#       Optionally allow a small kink (piecewise log-linear hazards) 
# with knot_age=75,
#       while taming curvature via lambda_kink.
#   (6) Winnow to 5-8 worlds by maximizing diversity (distance on log hazards)
#       subject to fit thresholds, and ensuring both system types represented.
# ------------------------------------------------------------------------------
# Major differences from the original:
# (1) parameters are a vector for easier optimization (this was a spec)
# (2) functions separate from their usage and updating (this was a spec)
# (3) ground world tweaked a bit
# (4) world creation is 2-steps now. We first perturb parameters in different 
#     directions using information from the svd of a jacobian of parameters on
#     lx and prev... it's knarly I admit. This replaces parameter forcing in the 
#     first version. THEN we optimize ALL parameters, but penalizing distance 
#     from the original perturbation.
# (5) then we keep worlds based on goodness of fit to prev and lx AND 
#     diversity in resulting transitions.
# (6) finally, rather than log linear hazards, I went with piecewise log linear
#     just so that the no-returns worlds would fit prevalence a bit better.
# We end up with 8 total worlds. Which could be called Octavio's Octet...
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("R/simulation_functions.R")

# ---------
# Top-level settings
# ---------

age      <- 50:100
age_int  <- 1
knot_age <- 75 # for piecewise log linear hazards
init     <- c(H=1, U=0) # just keep this simple

lambda      <- 0.3   # ridge toward each candidate (identity preservation)
                     # smaller = better fit to lx and prev; larger (closer to 1)
                     # respects the null direction parameters more.
lambda_kink <- 1.0   # tame curvature: penalize (s2 - s1)^2; set 0 for fully free kink
step_grid   <- c(-0.2, 0, 0.2) # how far to push parameters in null directions

# Fit thresholds used for winnowing down worlds 
# (otherwise we end up with 18 of them)
rmse_prev_max <- 4e-3    # 
rmse_lx_max   <- 1.3e-3

# How many final worlds?
K_final <- 8
K_noreturn_min <- 2

# ---------
# (1) Ground world (start linear, but we convert to piecewise with s2==s1)
# We can choose a different ground world of course. I got these values by
# eyeballing. I wanted prevalence to not end near 1. 

ground_pars_vec <- c(
  i_hu = -5.5, i_hd = -8.5, i_uh = -1.5, i_ud = -4.5,
  s_hu =  0.06, s_hd =  0.10, s_uh = -0.07, s_ud =  0.07
)

# Convert to piecewise params (linear baseline: s2 == s1)
ground_pars_pw <- as_piecewise_pars(ground_pars_vec)

# Ground Sullivan targets
ground_mslt <- run_mslt(ground_pars_pw, 
                        age = age, 
                        init = init, 
                        age_int = age_int, 
                        knot_age = knot_age)
ground_summary <- 
  ground_mslt |> 
  select(age, lx, prevalence)

# Quick check plot
ground_summary |>
  pivot_longer(c(lx, prevalence), 
               names_to = "Sullivan", 
               values_to = "value") |>
  ggplot(aes(age, value, color = Sullivan)) + 
  geom_line() + 
  labs(title = "Ground Sullivan outputs")


# ---------
# Helpers to generate + polish candidate worlds
# ---------

make_candidates <- function(base_pars, free_names, V2, grid_steps) {
  grid <- expand.grid(h1 = grid_steps, h2 = grid_steps)
  list(
    grid = grid,
    theta0 = purrr::pmap(grid, function(h1, h2) {
      apply_directions(base_pars = base_pars, V = V2, steps = c(h1, h2), free_names = free_names)
    })
  )
}

polish_many <- function(theta0_list,
                        ground_summary,
                        free_polish,
                        fixed = NULL,
                        lambda = 0.3,
                        lambda_kink = 0,
                        slope_floor = 1e-6) {
  
  purrr::map(theta0_list, ~ polish_candidate_ridge(
    theta_start = .x,
    theta_target = .x,
    ground_summary = ground_summary,
    free_polish = free_polish,
    fixed = fixed,
    age = age,
    init = init,
    age_int = age_int,
    knot_age = knot_age,
    lambda = lambda,
    lambda_kink = lambda_kink,
    slope_floor = slope_floor,
    maxit = 800
  ))
}

stack_mslt <- function(theta_list, meta_df = NULL, system = NA_character_) {
  purrr::imap_dfr(theta_list, function(th, i) {
    out <- run_mslt(th, age=age, init=init, age_int=age_int, knot_age=knot_age) |>
      mutate(pert = i, system = system)
    if (!is.null(meta_df)) {
      out <- out |>
        mutate(h1 = meta_df$h1[i], h2 = meta_df$h2[i])
    }
    out
  })
}

# ---------
# (2) RETURNS system: find null directions and generate candidates
# ---------

free_returns <- ms_par_names(piecewise = TRUE)

nd_ret <- find_null_directions(
  base_pars  = ground_pars_pw,
  free_names = free_returns,
  age = age,
  init = init,
  age_int = age_int,
  knot_age = knot_age,
  outputs = c("lx","prevalence"),
  tol = 1e-3,
  eps = 1e-4
)

V2_ret <- nd_ret$directions[, 1:2, drop = FALSE]

cand_ret <- make_candidates(ground_pars_pw, free_returns, V2_ret, step_grid)

# Ridge polish (all params) with kink penalty to keep curvature minor
theta_ret <- polish_many(
  theta0_list = cand_ret$theta0,
  ground_summary = ground_summary,
  free_polish = free_returns,
  fixed = NULL,
  lambda = lambda,
  lambda_kink = lambda_kink,
  slope_floor = 1e-6
)

# ---------
# (3) NO-RETURNS system: force UH ~ 0 and repeat
# ---------

uh_log <- log(1e-12)
fixed_noreturn <- c(i_uh = uh_log, s1_uh = 0, s2_uh = 0)

free_noreturn <- setdiff(free_returns, names(fixed_noreturn))

# base pars for no-return: start from ground but enforce UH fixed
base_nr <- ground_pars_pw
base_nr[names(fixed_noreturn)] <- fixed_noreturn

nd_nr <- find_null_directions(
  base_pars  = base_nr,
  free_names = free_noreturn,
  age = age,
  init = init,
  age_int = age_int,
  knot_age = knot_age,
  outputs = c("lx","prevalence"),
  tol = 1e-3,
  eps = 1e-4
)

V2_nr <- nd_nr$directions[, 1:2, drop = FALSE]

cand_nr <- make_candidates(base_nr, free_noreturn, V2_nr, step_grid)

# Enforce UH fixed in each candidate and ridge polish remaining parameters
cand_nr_theta0 <- purrr::map(cand_nr$theta0, function(th) {
  th[names(fixed_noreturn)] <- fixed_noreturn
  th
})

theta_nr <- polish_many(
  theta0_list = cand_nr_theta0,
  ground_summary = ground_summary,
  free_polish = free_noreturn,
  fixed = fixed_noreturn,
  lambda = lambda,
  lambda_kink = lambda_kink,
  slope_floor = 1e-6
)

# ---------
# (4) Combine, score fit, then select 5-8 diverse worlds

theta_all <- c(theta_ret, theta_nr)
system_all <- c(rep("returns", length(theta_ret)), rep("noreturn", length(theta_nr)))

fit_tab <- purrr::imap_dfr(theta_all, function(th, i) {
  score_fit(th, ground_summary, age, init=init, age_int=age_int, knot_age=knot_age) |>
    mutate(pert = i, system = system_all[i])
})

# had to eyeball limits. No returns-systems are harder to match
keep <- fit_tab |>
  filter(rmse_lx <= rmse_lx_max, rmse_prev <= rmse_prev_max)


keep_ids    <- keep$pert
theta_keep  <- theta_all[keep_ids]
system_keep <- system_all[keep_ids]

# Diversity distance on stacked log hazards
H <- do.call(rbind, lapply(theta_keep, haz_matrix, age=age, age_int=age_int, knot_age=knot_age))
D <- as.matrix(dist(H))

# Seed selection: one noreturn + one returns if available
idx_nr <- which(system_keep == "noreturn")
idx_rt <- which(system_keep == "returns")

seed <- integer(0)
if (length(idx_nr) > 0) seed <- c(seed, idx_nr[1])
if (length(idx_rt) > 0) seed <- c(seed, idx_rt[1])
if (length(seed) == 0) seed <- 1

# Greedy farthest-point selection with an initial set
sel <- unique(seed)
while (length(sel) < min(K_final, length(theta_keep))) {
  mind <- apply(D[, sel, drop=FALSE], 1, min)
  mind[sel] <- -Inf
  sel <- c(sel, which.max(mind))
}

# make sure at least 'K_noreturn_min' worlds are kept. We need
# this cuz their fits are a bit worse
if (K_noreturn_min > 0) {
  sel_sys <- system_keep[sel]
  need <- K_noreturn_min - sum(sel_sys == "noreturn")
  if (need > 0 && length(idx_nr) > 0) {
    cand_add <- setdiff(idx_nr, sel)
    if (length(cand_add) > 0) {
      # add the farthest noreturn candidates first
      order_add <- cand_add[order(apply(D[cand_add, sel, drop=FALSE], 1, min), decreasing = TRUE)]
      sel <- c(sel, head(order_add, need))
      sel <- sel[seq_len(min(length(sel), K_final))]
    }
  }
}

chosen_ids    <- keep_ids[sel]
theta_chosen  <- theta_all[chosen_ids]
system_chosen <- system_all[chosen_ids]

# here they are
print(fit_tab |> 
        filter(pert %in% chosen_ids) |> 
        arrange(system, rmse_prev))

# ---------
# (5) Stack outputs and plot
# ---------

lt_chosen <- purrr::imap_dfr(theta_chosen, function(th, i) {
  run_mslt(th, age=age, init=init, age_int=age_int, knot_age=knot_age) |>
    mutate(world = i, system = system_chosen[i])
})

# Sullivan outputs
lt_chosen |>
  ggplot(aes(age, lx, group=world, color=system)) +
  geom_line(alpha=.7) +
  geom_line(data=ground_summary, aes(age, lx), inherit.aes = FALSE, linewidth=1.1) +
  labs(title = "lx: chosen worlds vs ground",
       subtitle = "It's easy to match lx") +
  theme_minimal()

lt_chosen |>
  ggplot(aes(age, prevalence, group=world, color=system)) +
  geom_line(alpha=.7) +
  geom_line(data=ground_summary, aes(age, prevalence), inherit.aes = FALSE, linewidth=1.1) +
  labs(title = "Prevalence: chosen worlds vs ground",
       subtitle = "These variants are all within likely confidence bounds of survey-based prevalence") +
  theme_minimal()

# Hazards (log scale)
haz_chosen <- purrr::imap_dfr(theta_chosen, function(th, i) {
  expand_hazards(th, age=age, age_int=age_int, knot_age=knot_age) |>
    mutate(world=i, system=system_chosen[i])
})

haz_chosen |>
  ggplot(aes(age, hazard, group=world, color=system)) +
  geom_line(alpha=.7) +
  scale_y_log10() +
  facet_wrap(~trans, scales="free_y") +
  labs(title = "Hazards by chosen worlds (kink allowed, but kept small)") +
  theme_minimal()

haz_chosen |> 
  filter(trans %in% c("hd","ud")) |> 
  pivot_wider(names_from=trans,values_from=hazard) |> 
  mutate(Ra = ud / hd) |> 
  ggplot(aes(x=age,y=Ra,color = as.factor(world))) +
  geom_line() +
  ylim(1,80) +
  labs(title = "Risk ratios for different worlds")  +
  theme_minimal()

lt_chosen |> 
  select(world,age, HU,HD,UH,UD) |> 
  pivot_longer(-c(age,world), names_to = "trans",values_to = "p") |> 
  ggplot(aes(x=age,y=p,color=as.factor(world))) +
  geom_line() +
  facet_wrap(~trans, scales = "free_y") + 
  theme_minimal() +
  labs(title = "Variety in transitions over worlds")

# save world parameters (minimal representation)

bind_rows(theta_chosen) |> 
  mutate(system = system_chosen,
         world = 1:n(), .before=1) |> 
  write_csv("data/worlds_octet.csv")
