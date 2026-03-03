# ------------------------------------------------------------------------------
#
# Goal:
#   Build a small set (5-8) of multistate "worlds" that match the same Sullivan
#   outputs (lx, prevalence), including ~2-3 "no-return" worlds (UH ~ 0).
#
# Method (current):
#   (1) Define a ground world and compute its Sullivan targets (lx, prevalence).
#   (2) Find near-null directions (SVD of Jacobian of (lx, prevalence)).
#   (3) Step along (dir1, dir2) to create candidate parameter sets.
#   (4) Polish candidates with ridge-regularized bounded optimization (optional 2nd pass).
#   (5) From the *first-pass* polished returns worlds, extract Rx(age)=ud/hd curves.
#   (6) For each Rx curve, derive an *exact* RETURNS hazard schedule that matches
#       ground lx & prevalence by solving for (hu,uh) per age interval.
#   (7) Create 2 no-return worlds by deriving hazards directly from hand-picked Rx curves.
#   (8) Export hazards + diagnostics plots.
# ------------------------------------------------------------------------------

library(tidyverse)
source("R/simulation_functions3.R")

# ---------
# Top-level settings
# ---------

age     <- 50:100
age_int <- 1
init    <- c(H = 1, U = 0)

# Hazard shape used ONLY for the "candidate world" polishing stage (returns worlds)
pivot_age   <- 70     # location of slope transition (shared across transitions)
shape_width <- 6      # smoothness (years; used for softkink)
haz_shape   <- "softkink"  # "loglinear", "piecewise", "softkink"

# Regularization / fit control (returns polishing)
lambda1      <- 0.3    # ridge toward each candidate (identity preservation)
lambda2      <- 0.01   # optional second pass (smaller => tighter fit, but can shrink diversity)
do_pass2     <- F   # set FALSE to skip pass-2 entirely

lambda_kink1 <- 0.5
lambda_kink2 <- 0
kink_trans1  <- c("hd","ud")
kink_trans2  <- c("hd","ud")

w_lx  <- 0.5
w_prev <- 0.5
age_weight_lx   <- "none"
age_weight_prev <- "none"
prev_band       <- c(52, 80)
prev_band_mult  <- 1

maxit1 <- 800
maxit2 <- 1500
factr1 <- 1e6
factr2 <- 1e4

# Candidate stepping along null directions
step_grid <- c(-0.2, 0, 0.2)

# Fit thresholds for screening (applied to *polished* theta, before diversity selection)
rmse_prev_max <- 4e-3
rmse_lx_max   <- 1.3e-3

# Final set size
K_final    <- 8
K_noreturn <- 2

# Which world id is reserved for the ground world (keeps your tiny sanity test stable)
GROUND_WORLD_ID <- 5

# ---------
# (1) Ground world (log-linear baseline, represented as piecewise with s2==s1)
# ---------
ground_pars_vec <- c(
  i_hu = -5.5, i_hd = -8.5, i_uh = -1.5, i_ud = -4.5,
  s_hu =  0.06, s_hd =  0.10, s_uh = -0.07, s_ud =  0.07
)

ground_pars_pw <- as_piecewise_pars(ground_pars_vec)

ground_mslt <- run_mslt(
  ground_pars_pw,
  age = age, init = init, age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
)

ground_summary <- ground_mslt |> select(age, lx, prevalence)

ground_summary |>
  pivot_longer(c(lx, prevalence), names_to = "Sullivan", values_to = "value") |>
  ggplot(aes(age, value, color = Sullivan)) +
  geom_line() +
  labs(title = "Ground Sullivan outputs")

# ---------
# (2) RETURNS system: null directions + candidates + polish
# ---------
free_returns <- ms_par_names(piecewise = TRUE)

nd_ret <- find_null_directions(
  base_pars  = ground_pars_pw,
  free_names = free_returns,
  age = age, init = init, age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width,
  outputs = c("lx","prevalence"),
  tol = 1e-3, eps = 1e-4
)

V2_ret <- nd_ret$directions[, 1:2, drop = FALSE]

cand_ret <- make_candidates(
  base_pars  = ground_pars_pw,
  free_names = free_returns,
  V2 = V2_ret,
  grid_steps = step_grid
)

# Pass 1 (always): used for Rx extraction (keeps diversity)
theta_ret_1 <- polish_many(
  theta0_list = cand_ret$theta0,
  ground_summary = ground_summary,
  free_polish = free_returns,
  age = age, init = init, age_int = age_int,
  pivot_age = pivot_age, shape_width = shape_width,
  lambda = lambda1,
  lambda2 = NULL, # <- IMPORTANT: first pass only
  lambda_kink1 = lambda_kink1, kink_trans1 = kink_trans1,
  lambda_kink2 = lambda_kink2, kink_trans2 = kink_trans2,
  maxit1 = maxit1, maxit2 = maxit2,
  factr1 = factr1, factr2 = factr2,
  w_lx = w_lx, w_prev = w_prev,
  age_weight_lx = age_weight_lx, age_weight_prev = age_weight_prev,
  prev_band = prev_band, prev_band_mult = prev_band_mult,
  n_cores = 4
)

# Optional pass 2: tighter fit but can shrink diversity
theta_ret <- if (isTRUE(do_pass2)) {
  polish_many(
    theta0_list = theta_ret_1,          # warm start from pass 1
    ground_summary = ground_summary,
    free_polish = free_returns,
    age = age, init = init, age_int = age_int,
    pivot_age = pivot_age, shape_width = shape_width,
    lambda = lambda2,                   # second-stage lambda
    lambda2 = NULL,                     # no 3rd stage
    lambda_kink1 = lambda_kink2, kink_trans1 = kink_trans2,
    lambda_kink2 = 0, kink_trans2 = kink_trans2,
    maxit1 = maxit2, maxit2 = 0,
    factr1 = factr2, factr2 = factr2,
    w_lx = w_lx, w_prev = w_prev,
    age_weight_lx = age_weight_lx, age_weight_prev = age_weight_prev,
    prev_band = prev_band, prev_band_mult = prev_band_mult,
    n_cores = 4
  )
} else {
  theta_ret_1
}

# ---------
# (3) Select RETURNS worlds (based on fit + diversity), then extract Rx from *pass 1*
# ---------
fit_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  score_fit(
    th, ground_summary,
    age = age, init = init, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |> mutate(pert = i)
})

keep_ret <- fit_ret |>
  filter(rmse_lx <= rmse_lx_max, rmse_prev <= rmse_prev_max)

theta_pool <- theta_ret[keep_ret$pert]
theta_pool_1 <- theta_ret_1[keep_ret$pert]   # <-- SAME indices, but pass-1 versions

if (length(theta_pool) < (K_final - 1)) {
  warning("Too few returns worlds passed thresholds; using best-fitting ones instead.", call. = FALSE)
  ord <- order(fit_ret$rmse_prev + fit_ret$rmse_lx)
  take <- ord[seq_len(min(length(theta_ret), K_final - 1))]
  theta_pool <- theta_ret[take]
  theta_pool_1 <- theta_ret_1[take]
}

# Diversity selection among the pool (excluding ground, which we force)
H_pool <- do.call(rbind, lapply(theta_pool_1, haz_feat_from_theta))
D_pool <- as.matrix(dist(H_pool))
sel_pool <- select_farthest(D_pool, k = K_final - 1)

theta_sel   <- theta_pool[sel_pool]
theta_sel_1 <- theta_pool_1[sel_pool]

# Place the ground world at fixed ID
theta_chosen_1 <- vector("list", K_final)   # pass-1 thetas (for Rx extraction)
theta_chosen   <- vector("list", K_final)   # polished thetas (for reference only)

theta_chosen_1[[GROUND_WORLD_ID]] <- ground_pars_pw
theta_chosen[[GROUND_WORLD_ID]]   <- ground_pars_pw

fill_slots <- setdiff(seq_len(K_final), GROUND_WORLD_ID)
theta_chosen_1[fill_slots] <- theta_sel_1[seq_len(min(length(theta_sel_1), length(fill_slots)))]
theta_chosen[fill_slots]   <- theta_sel[seq_len(min(length(theta_sel), length(fill_slots)))]

# Hazards for chosen returns worlds (from PASS 1; these define Rx patterns)
haz_ret_pass1 <- purrr::imap_dfr(theta_chosen_1, function(th, i) {
  expand_hazards(
    pars = th, age = age, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(world = i, system = "returns_pass1")
})

Rx_plot_df <- haz_ret_pass1 |>
  filter(trans %in% c("hd","ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd) |>
  select(age, world, Rx)

ggplot(Rx_plot_df, aes(x = age, y = Rx, color = as.factor(world))) +
  geom_line() +
  labs(color = "returns world", y = "Rx = ud/hd",
       title = "Rx curves from chosen returns worlds (PASS 1; pick low/high for no-returns)")

# ---------
# (4) Build EXACT RETURNS hazards from Rx (matches lx & prevalence by construction)
# ---------

haz_returns_exact <- purrr::map_dfr(split(Rx_plot_df, Rx_plot_df$world), function(dfw) {
  w <- unique(dfw$world)
  Rx_vec <- dfw |> dplyr::arrange(age) |> dplyr::pull(Rx)
  
  derive_returns_hazards_from_Rx(
    age  = ground_summary$age,
    lx   = ground_summary$lx,
    prev = ground_summary$prevalence,
    Rx   = Rx_vec,
    age_int = age_int,
    bounds = c(1e-12, 2),
    maxit = 500,
    reltol = 1e-14,
    smooth_w = 0.5,     # try 0, 0.2, 0.5, 1
    bound_w  = 1e-3,    # try 0, 1e-4, 1e-3
    verbose = FALSE
  ) |>
    dplyr::mutate(world = w, system = "returns")
})
haz_returns_exact |> ggplot(aes(x=age,y=hazard, color = as.factor(world)))+ 
  geom_line() + 
  facet_wrap(~trans) +
  scale_y_log10()
lt_exact <- haz_returns_exact |>
  group_split(world) |>
  map_dfr(function(hz) {
    w <- unique(hz$world)
    sys <- unique(hz$system)
    P <- haz_to_probs(hz |> select(age, trans, hazard), age = age, age_int = age_int)
    calculate_lt(P, init = init) |> mutate(world = w, system = sys)
  })
lt_exact |> 
  select(world, age, lx) |> 
  pivot_wider(names_from = world, values_from = lx) |> 
  mutate(base_pattern = `5`) |> 
  pivot_longer(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`), names_to = "world", values_to = "lx") |> 
  mutate(residual = base_pattern - lx) |> 
  ggplot(aes(x = age,y=residual, color = world)) +
  geom_line() +
  labs(title = "lx residuals")

lt_chosen |> 
  select(world, age, prevalence) |> 
  pivot_wider(names_from = world, values_from = prevalence) |> 
  mutate(base_pattern = `5`) |> 
  pivot_longer(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`), names_to = "world", values_to = "prevalence") |> 
  mutate(residual = base_pattern - prevalence) |> 
  ggplot(aes(x = age,y=residual, color = world)) +
  geom_line()+
  labs(title = "prev residuals")

# ---------
# (5) NO-RETURNS worlds (hand-pick two Rx sources, derive hazards, swap in)
# ---------
Rx_source_worlds <- c(3, 4)  # <-- YOU eyeball these from the Rx plot (low/high)
nr_world_ids     <- c(1, 6)  # which world IDs to replace with no-returns

stopifnot(length(Rx_source_worlds) == length(nr_world_ids))
stopifnot(length(Rx_source_worlds) == K_noreturn)

haz_nr_list <- purrr::map2(Rx_source_worlds, nr_world_ids, function(src_w, out_w) {
  Rx_vec <- Rx_plot_df |>
    filter(world == src_w) |>
    arrange(age) |>
    pull(Rx)
  
  derive_noreturns_hazards(
    age = ground_summary$age,
    lx  = ground_summary$lx,
    prev = ground_summary$prevalence,
    Rx  = Rx_vec,
    age_int = age_int,
    verbose = FALSE
  ) |>
    mutate(world = out_w, system = "noreturn")
})

haz_nr <- bind_rows(haz_nr_list)

# Replace selected world IDs with their no-return versions
haz_chosen <- haz_returns_exact |>
  filter(!(world %in% nr_world_ids)) |>
  bind_rows(haz_nr) |>
  arrange(world, trans, age)

# ---------
# (6) Export + quick plots
# ---------
dir.create("data", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(haz_chosen, "data/haz_octet.csv.gz")

haz_chosen |>
  filter(hazard > 1e-10) |>
  ggplot(aes(x = age, y = hazard, color = as.factor(world), linetype = system)) +
  geom_line(alpha = .85) +
  facet_wrap(~trans, scales = "free_y") +
  scale_y_log10() +
  labs(color = "world", title = "Hazards for chosen worlds (returns + no-returns)") +
  theme_minimal()

# ---------
# (7) Tiny sanity test: world 5 is the ground world
# ---------
lt_world5 <- calculate_lt(
  haz_to_probs(
    haz_chosen |> filter(world == GROUND_WORLD_ID) |> select(age, trans, hazard),
    age = age, age_int = age_int
  ),
  init = init
)

print(max(abs(lt_world5$lx - ground_summary$lx)))
print(max(abs(lt_world5$prevalence - ground_summary$prevalence)))

# Sullivan outputs for all chosen worlds
lt_chosen <- haz_chosen |>
  group_split(world) |>
  map_dfr(function(hz) {
    w <- unique(hz$world)
    sys <- unique(hz$system)
    P <- haz_to_probs(hz |> select(age, trans, hazard), age = age, age_int = age_int)
    calculate_lt(P, init = init) |> mutate(world = w, system = sys)
  })

# Residual plots (vs ground world)
lt_chosen |>
  select(world, age, lx) |>
  pivot_wider(names_from = world, values_from = lx) |>
  mutate(base = .data[[as.character(GROUND_WORLD_ID)]]) |>
  pivot_longer(-c(age, base), names_to = "world", values_to = "lx") |>
  mutate(residual = base - lx) |>
  ggplot(aes(x = age, y = residual, color = world)) +
  geom_line() +
  labs(title = "lx residuals (vs ground)")

lt_chosen |>
  select(world, age, prevalence) |>
  pivot_wider(names_from = world, values_from = prevalence) |>
  mutate(base = .data[[as.character(GROUND_WORLD_ID)]]) |>
  pivot_longer(-c(age, base), names_to = "world", values_to = "prevalence") |>
  mutate(residual = base - prevalence) |>
  ggplot(aes(x = age, y = residual, color = world)) +
  geom_line() +
  labs(title = "prevalence residuals (vs ground)")

hle_sullivan <- lt_chosen |>
  group_by(world) |>
  summarize(hle_sullivan = sum((1 - prevalence) * lx), .groups = "drop") |>
  pull(hle_sullivan)

round(hle_sullivan - hle_sullivan[GROUND_WORLD_ID], 8)