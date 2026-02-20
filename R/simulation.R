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
    theta0 = pmap(grid, function(h1, h2) {
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
  
  map(theta0_list, ~ polish_candidate_ridge(
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
  imap_dfr(theta_list, function(th, i) {
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
# (3) Select diverse RETURNS worlds (null directions + ridge polishing already done)
# ---------

fit_tab_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  score_fit(th, ground_summary, age, init = init, age_int = age_int, knot_age = knot_age) |>
    mutate(pert = i, system = "returns")
})

keep_ret <- fit_tab_ret |>
  filter(rmse_lx <= rmse_lx_max, rmse_prev <= rmse_prev_max)

keep_ids_ret <- keep_ret$pert
theta_keep   <- theta_ret[keep_ids_ret]

# Diversity distance on stacked log hazards (returns only)
H_ret <- do.call(rbind, lapply(theta_keep, haz_matrix, age = age, age_int = age_int, knot_age = knot_age))
D_ret <- as.matrix(dist(H_ret))

# Greedy farthest-point selection (seed with the most central keep, i.e. medoid)
seed <- if (nrow(D_ret) > 0) which.min(rowSums(D_ret)) else 1

sel <- seed
while (length(sel) < min(K_final, length(theta_keep))) {
  mind <- apply(D_ret[, sel, drop = FALSE], 1, min)
  mind[sel] <- -Inf
  sel <- c(sel, which.max(mind))
}

theta_chosen <- theta_keep[sel]

# Hazards for chosen RETURNS worlds (this is the primary representation going forward)
haz_chosen <-
  imap_dfr(theta_chosen, function(th, i) {
    expand_hazards(th, age = age, age_int = age_int, knot_age = knot_age) |>
      mutate(world = i, system = "returns")
  }) |>
  arrange(world, trans, age)

# Quick summary of fits for the chosen RETURNS worlds
print(
  keep_ret |>
    mutate(world = match(pert, keep_ids_ret)) |>
    filter(pert %in% keep_ids_ret[sel]) |>
    arrange(rmse_prev)
)

# ---------
# (4) NO-RETURNS worlds: derive hazards directly from hand-picked Rx patterns
#     (no null-direction search / no optimization)
# ---------

# Plot Rx(age)=ud/hd for the chosen RETURNS worlds, then set the IDs below by eyeballing
Rx_plot_df <-
  haz_chosen |>
  filter(trans %in% c("hd", "ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd)

ggplot(Rx_plot_df, ggplot2::aes(x = age, y = Rx, color = as.factor(world))) +
  geom_line(linewidth = 1.0) +
  labs(color = "returns world", y = "Rx = ud/hd",
       title = "Choose low/high Rx curves from RETURNS worlds") +
  theme_minimal()

# After eyeballing the plot, set:
#   Rx_source_worlds = which RETURNS worlds provide Rx patterns (e.g., low/high)
#   nr_world_ids     = which world IDs in the final set you want to become NO-RETURNS
Rx_source_worlds <- c(3, 4)  # 
nr_world_ids     <- c(1, 6)  # 

stopifnot(length(Rx_source_worlds) == length(nr_world_ids))

haz_nr_list <- map2(Rx_source_worlds, nr_world_ids, function(src_w, out_w) {
  
  Rx_vec <- Rx_plot_df |>
    filter(world == src_w) |>
    arrange(age) |>
    pull(Rx)
  
  derive_noreturns_hazards(
    age  = ground_summary$age,
    lx   = ground_summary$lx,
    prev = ground_summary$prevalence,
    Rx   = Rx_vec,
    age_int = age_int,
    verbose = FALSE
  ) |>
    mutate(world = out_w, system = "noreturn")
})

haz_nr <- bind_rows(haz_nr_list)

# Replace those worlds with their no-return counterparts
haz_chosen <-
  haz_chosen |>
  filter(!(system == "returns" & world %in% nr_world_ids)) |>
  bind_rows(haz_nr) |>
  arrange(world, trans, age)

# ---------
# (5) Sullivan outputs + plots + exports (hazards are primary)
# ---------

# Build lifetable outputs from hazards for ALL chosen worlds (returns and no-returns)
lt_chosen <-
  haz_chosen |>
  group_split(world) |>
  map_dfr(function(hz) {
    w <- unique(hz$world)
    sys <- unique(hz$system)
    P <- haz_to_probs(hz |> select(age, trans, hazard), age = age, age_int = age_int)
    calculate_lt(P, init = init) |>
      mutate(world = w, system = sys)
  })

# Sullivan outputs
lt_chosen |>
  ggplot(aes(age, lx, group = world, color = as.factor(world))) +
  geom_line(alpha = .75) +
  geom_line(data = ground_summary, ggplot2::aes(age, lx),
            inherit.aes = FALSE, linewidth = 1.1) +
  labs(title = "lx: chosen worlds vs ground") +
  theme_minimal()

lt_chosen |>
  ggplot(aes(age, prevalence, group = world, color = as.factor(world))) +
  geom_line(alpha = .75) +
  geom_line(data = ground_summary, ggplot2::aes(age, prevalence),
            inherit.aes = FALSE, linewidth = 1.1) +
  labs(title = "Prevalence: chosen worlds vs ground") +
  theme_minimal()

# Hazards (log scale)
haz_chosen |>
  filter(hazard > 1e-10) |> 
  ggplot(aes(age, hazard, color = as.factor(world), linetype = system)) +
  geom_line(alpha = .85) +
  scale_y_log10() +
  facet_wrap(~trans, scales = "free_y") +
  labs(title = "Hazards by chosen worlds (includes no-returns)") +
  theme_minimal()

# Risk ratios
haz_chosen |>
  filter(trans %in% c("hd", "ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd) |>
  ggplot(aes(age, Rx, color = as.factor(world), linetype = system)) +
  geom_line(linewidth = 1.1) +
  labs(title = "Rx = ud/hd by world",
       subtitle = "Note the no-returns RRs overlap other ones, as these were generative") +
  theme_minimal()

# Save worlds
write_csv(haz_chosen, "data/haz_octet.csv.gz")
