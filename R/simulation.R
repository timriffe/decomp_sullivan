# simulation.R
# Cleaned application script for multistate 'world' generation + exact snapping.
#
# Pipeline:
#   1) Define ground world, compute Sullivan targets (lx, prevalence).
#   2) Find near-null directions via SVD of Jacobian of (lx, prevalence).
#   3) Create candidate parameter sets stepping along first two near-null directions
#      (excludes the (0,0) step by construction).
#   4) Single-pass polish candidates with polish_many() (no second-stage refinement).
#   5) Extract Rx(age)=ud/hd curves from polished candidates.
#   6) For each Rx curve, compute exact RETURNS hazards by isocline snapping of (hu0,uh0).
#   7) For each Rx curve, compute deterministic NO-RETURNS hazards.
#   8) Bind hazards and compute lifetables via group_modify() (no purrr::map in this step).
#
# You should be able to run this script after sourcing simulation_functions.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/simulation_functions5.R")

# -----------------------------------------------------------------------------
# USER SETTINGS (keep these minimal)
# -----------------------------------------------------------------------------
age     <- 50:100
age_int <- 1
init    <- c(H = 1, U = 0)

# Hazard family for ground + candidates + polishing
haz_shape   <- "softkink"   # "loglinear", "piecewise", "softkink"
pivot_age   <- 70
shape_width <- 6

# Candidate stepping along null directions (NOTE: (0,0) excluded automatically)
step_grid <- c(-0.3, -.1, .1, 0.3)   
n_cores=parallel::detectCores()
# Single-pass polishing controls
lambda1 <- 0.3
maxit1  <- 800
factr1  <- 1e6
w_lx    <- 0.5
w_prev  <- 0.5

# Screening thresholds (optional; set to Inf to keep everything)
rmse_prev_max <- Inf
rmse_lx_max   <- Inf

# Exact snapping controls
bounds_scalar <- c(1e-12, 2)
turnover_K <- 1.2
snap_tol  <- 1e-20
line_tol  <- 1e-12
eps_log   <- 1e-2

# ---------------------------------------------------
# 1) Ground world
# ---------------------------------------------------
# This arbitrary parameterization creates the ground world.
# from here onward we try to match the sullivan stats
ground_pars_vec <- c(
  i_hu = -5.5, i_hd = -8.5, i_uh = -1.5, i_ud = -4.5,
  s_hu =  0.06, s_hd =  0.10, s_uh = -0.07, s_ud =  0.07
)
ground_pars <- as_piecewise_pars(ground_pars_vec)

ground_mslt <- run_mslt(
  ground_pars,
  age = age, init = init, age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
)

ground_summary <- ground_mslt |> select(age, lx, prevalence)

# Quick plot of targets
ground_summary |>
  pivot_longer(c(lx, prevalence), names_to = "measure", values_to = "value") |>
  ggplot(aes(age, value, color = measure)) +
  geom_line() +
  labs(title = "Ground Sullivan targets")

# ---------------------------------------
# 2) Null directions from Jacobian
# ---------------------------------------
free_names <- ms_par_names(piecewise = TRUE)

nd <- find_null_directions(
  base_pars  = ground_pars,
  free_names = free_names,
  age = age, init = init, age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width,
  outputs = c("lx","prevalence"),
  tol = 1e-3, eps = 1e-4
)

V2 <- nd$directions[, 1:2, drop = FALSE]

# -------------------------------------------------------
# 3) Candidate parameter sets along V2 (exclude (0,0))
# -------------------------------------------------------
cand <- make_candidates(
  base_pars  = ground_pars,
  free_names = free_names,
  V2 = V2,
  grid_steps = step_grid,
  exclude_zero = TRUE
)

# ---------------------------------------------
# 4) Single-pass polish
# ---------------------------------------------
theta_polished <- polish_many(
  theta0_list = cand$theta0,
  ground_summary = ground_summary,
  free_polish = free_names,
  age = age, init = init, age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width,
  lambda = lambda1,
  lambda2 = NULL,            # IMPORTANT: no 2nd stage
  lambda_kink1 = 0,
  kink_trans1 = character(0),
  maxit1 = maxit1,
  factr1 = factr1,
  w_lx = w_lx,
  w_prev = w_prev,
  age_weight_lx = "none",
  age_weight_prev = "none",
  prev_band = NULL,
  prev_band_mult = 1,
  n_cores = n_cores
)

fit_tbl <- purrr::imap_dfr(theta_polished, function(th, i){
  score_fit(
    th, ground_summary,
    age = age, init = init, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |> mutate(world = i)
})

keep_idx <- fit_tbl |>
  filter(.data$rmse_lx <= rmse_lx_max, .data$rmse_prev <= rmse_prev_max) |>
  pull(.data$world)

if(length(keep_idx) == 0){
  warning("No polished worlds passed thresholds; keeping all.", call. = FALSE)
  keep_idx <- seq_along(theta_polished)
}

theta_keep <- theta_polished[keep_idx]

# Hazards for kept returns worlds (pass-1 polished)
haz_ret_pass1 <- purrr::imap_dfr(theta_keep, function(th, j){
  expand_hazards(
    pars = th, age = age, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(world = keep_idx[j], system = "returns_pass1")
})

# -------------------------------
# 5) Extract Rx = ud/hd
# -------------------------------
Rx_plot_df <- haz_ret_pass1 |>
  filter(.data$trans %in% c("hd","ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd) |>
  select(age, world, Rx)

ggplot(Rx_plot_df, aes(age, Rx, color = as.factor(world))) +
  geom_line() +
  labs(color = "world", y = "Rx = ud/hd", title = "Rx curves from polished candidates")

# ------------------------------------------------
# 6) Exact RETURNS hazards by isocline snapping
# ------------------------------------------------
haz_returns_exact <- Rx_plot_df %>%
  group_by(.data$world) %>%
  group_modify(~{
    w <- unique(.y$world)
    Rx_vec <- .x %>% arrange(age) %>% pull(Rx)
    
    hz0 <- haz_ret_pass1 %>%
      filter(.data$world == w, .data$trans %in% c("hu","uh")) %>%
      arrange(.data$trans, .data$age)
    
    hu0 <- hz0 %>% filter(.data$trans == "hu") %>% pull(.data$hazard)
    uh0 <- hz0 %>% filter(.data$trans == "uh") %>% pull(.data$hazard)
    
    derive_returns_hazards_from_Rx(
      age  = ground_summary$age,
      lx   = ground_summary$lx,
      prev = ground_summary$prevalence,
      Rx   = Rx_vec,
      age_int = age_int,
      bounds = bounds_scalar,
      hu_bounds = NULL,
      uh_bounds = NULL,
      hu_0 = hu0,
      uh_0 = uh0,
      turnover_K = turnover_K,
      snap_tol = snap_tol,
      line_tol = line_tol,
      eps_log = eps_log,
      verbose = FALSE
    ) %>% mutate(system = "returns")
  }) %>%
  ungroup()

# ------------------------------
# 7) No-returns hazards
# ------------------------------
haz_noreturns <- Rx_plot_df %>%
  group_by(.data$world) %>%
  group_modify(~{
    w <- unique(.y$world)
    Rx_vec <- .x %>% arrange(age) %>% pull(Rx)
    
    derive_noreturns_hazards(
      age  = ground_summary$age,
      lx   = ground_summary$lx,
      prev = ground_summary$prevalence,
      Rx   = Rx_vec,
      age_int = age_int,
      verbose = FALSE
    ) %>% mutate(system = "noreturns")
  }) %>%
  ungroup()

haz_all <- bind_rows(haz_returns_exact, haz_noreturns)

# -----------------------------------
# 8) Lifetables from hazards 
# -----------------------------------
lt_all <- haz_all %>%
  group_by(.data$world, .data$system) %>%
  group_modify(~{
    P <- haz_to_probs(.x %>% select(age, trans, hazard), age = age, age_int = age_int)
    calculate_lt(P, init = init)
  }) %>%
  ungroup()


lt_all |> 
  ggplot(aes(x=age,y=prevalence,color=interaction(world,system))) +
  geom_line()

haz_all |> 
  filter(!(trans=="uh" & system== "noreturns")) |> 
  ggplot(aes(x=age,y=hazard,color=as.factor(world))) +
  geom_line()+
  scale_y_log10()+
  facet_grid(vars(trans),vars(system),scale="free_y")
  

lt_all |> 
  group_by(system,world) |> 
  summarize(hle = sum(Lx * (1-prevalence))) |> 
  View()
