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
#       Optionally allow a smooth slope change in log-hazards (softkink) with
#       pivot_age=75, while taming curvature via lambda_kink.
#   (6) Winnow to 5-8 worlds by maximizing diversity (distance on log hazards)
#       subject to fit thresholds.
#   (7) Create 2-3 no-return worlds by *deriving hazards directly* from hand-
#       picked Rx(age)=ud/hd curves (no SVD/optim for no-returns), and swap them
#       into the final set.
# ------------------------------------------------------------------------------
# We end up with 8 total worlds. Which could be called Octavio's Octet...
library(tidyverse)

source("R/simulation_functions2.R")

# ---------
# Top-level settings
# ---------

age      <- 50:100
age_int  <- 1
init <- c(H=1, U=0)
# Hazard shape:

pivot_age   <- 70   # location of slope transition (shared across transitions)
shape_width <- 6    # smoothness of slope transition (years; used for softkink)
haz_shape   <- "softkink" # "loglinear", "piecewise"
lambda1      <- 0.3   # ridge toward each candidate (identity preservation)
lambda2      <- 0.01
lambda_kink1 <- .5   # tame curvature: penalize (s2 - s1)^2; set 0 for fully free kink
lambda_kink2 <- 0
kink_trans1 = c("hd","ud")
kink_trans2 = c("hd","ud")
w_lx = 0.5
w_prev = .5
age_weight_lx = "none"
age_weight_prev = "none"
prev_band = c(52, 80)
prev_band_mult = 1
maxit1 = 800
maxit2 = 1500

factr1 = 1e6
factr2 = 1e4
w_int = 10 
integrand = "health"
age_weight_int = "none"
# Candidate stepping along null directions
step_grid   <- c(-0.2, 0, 0.2)

# Fit thresholds used for winnowing down worlds
rmse_prev_max <- 4e-3
rmse_lx_max   <- 1.3e-3

# How many final worlds?
K_final <- 8

# How many no-return worlds to swap in (hand-picked)
K_noreturn <- 2

# ---------
# (1) Ground world (start linear, then convert to piecewise with s2==s1)
# ---------

ground_pars_vec <- c(
  i_hu = -5.5, i_hd = -8.5, i_uh = -1.5, i_ud = -4.5,
  s_hu =  0.06, s_hd =  0.10, s_uh = -0.07, s_ud =  0.07
)

# Convert to piecewise params (linear baseline: s2 == s1)
ground_pars_pw <- as_piecewise_pars(ground_pars_vec)

# Ground Sullivan targets
ground_mslt <- run_mslt(
  ground_pars_pw,
  age = age,
  init = init,
  age_int = age_int,
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
)

ground_summary <- ground_mslt |> select(age, lx, prevalence)

# Quick check plot
ground_summary |>
  pivot_longer(c(lx, prevalence), names_to = "Sullivan", values_to = "value") |>
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
      apply_directions(
        base_pars = base_pars,
        V = V2,
        steps = c(h1, h2),
        free_names = free_names
      )
    })
  )
}

polish_many <- function(theta0_list,
                        ground_summary,
                        free_polish,
                        age,
                        init = c(H=1,U=0),
                        age_int = 1,
                        pivot_age = 75,
                        shape_width = 5,
                        knot_age = NULL,
                        fixed = NULL,
                        lambda = 0.3,
                        lambda2 = NULL,
                        slope_floor = 1e-6,
                        maxit1 = 800,
                        maxit2 = 400,
                        # kink pass-specific
                        lambda_kink1 = 0,
                        kink_trans1 = c("hu","hd","uh","ud"),
                        lambda_kink2 = lambda_kink1,
                        kink_trans2 = kink_trans1,
                        # fit weights
                        w_lx = 1,
                        w_prev = 1,
                        age_weight_lx = c("none","lx","Lx"),
                        age_weight_prev = c("none","lx","Lx"),
                        prev_band = c(52, 70),
                        prev_band_mult = 1,
                        # integrand
                        w_int = 0,
                        integrand = c("health","unhealthy"),
                        age_weight_int = c("none","lx","Lx"),
                        # convergence
                        factr1 = 1e6,
                        factr2 = 1e6,
                        # parallel
                        n_cores = 1) {
  
  age_weight_lx <- match.arg(age_weight_lx)
  age_weight_prev <- match.arg(age_weight_prev)
  integrand <- match.arg(integrand)
  age_weight_int <- match.arg(age_weight_int)
  
  one <- function(th0) {
    polish_candidate_ridge(
      theta_start = th0,
      theta_target = th0,
      ground_summary = ground_summary,
      free_polish = free_polish,
      age = age,
      init = init,
      age_int = age_int,
      haz_shape = haz_shape,  # matches your script style
      pivot_age = pivot_age,
      shape_width = shape_width,
      knot_age = knot_age,
      fixed = fixed,
      lambda = lambda,
      lambda2 = lambda2,
      slope_floor = slope_floor,
      maxit = maxit1,
      maxit2 = maxit2,
      lambda_kink1 = lambda_kink1,
      kink_trans1 = kink_trans1,
      lambda_kink2 = lambda_kink2,
      kink_trans2 = kink_trans2,
      w_lx = w_lx,
      w_prev = w_prev,
      age_weight_lx = age_weight_lx,
      age_weight_prev = age_weight_prev,
      prev_band = prev_band,
      prev_band_mult = prev_band_mult,
      w_int = w_int,
      integrand = integrand,
      age_weight_int = age_weight_int,
      factr1 = factr1,
      factr2 = factr2
    )
  }
  
  if (is.null(n_cores) || !is.finite(n_cores) || n_cores <= 1) {
    return(purrr::map(theta0_list, one))
  }
  
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Parallel polish_many() requires packages 'future' and 'future.apply'. Install them or set n_cores = 1.",
         call. = FALSE)
  }
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = as.integer(n_cores))
  
  future.apply::future_lapply(theta0_list, one, future.seed = TRUE)
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
  haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width,
  outputs = c("lx","prevalence"),
  tol = 1e-3,
  eps = 1e-4
)

V2_ret <- nd_ret$directions[, 1:2, drop = FALSE]
cand_ret <- make_candidates(ground_pars_pw, free_returns, V2_ret, step_grid)

# Polish (all params) with kink penalty to keep curvature minor
# stopifnot(is.list(cand_ret$theta0))
# stopifnot(is.numeric(cand_ret$theta0[[1]]))
# stopifnot(!is.null(names(cand_ret$theta0[[1]])))
# 
# setdiff(free_returns, names(cand_ret$theta0[[1]]))

theta_ret <- polish_many(
  theta0_list = cand_ret$theta0,
  ground_summary = ground_summary,
  free_polish = free_returns,
  age = age,
  init = init,
  age_int = age_int,
  pivot_age = pivot_age,
  lambda = lambda1,
  lambda2 = lambda2,
  lambda_kink1 = lambda_kink1, 
  kink_trans1 = kink_trans1,
  lambda_kink2 = lambda_kink2, 
  kink_trans2 = kink_trans2, 
  maxit1 = maxit1,
  maxit2 = maxit2,
  factr1 = factr1,
  factr2 = factr2,
  w_lx = w_lx,
  w_int = w_int, 
  integrand = integrand,
  w_prev = w_prev,
  age_weight_lx = age_weight_lx,
  age_weight_prev = age_weight_prev,
  prev_band = prev_band,
  prev_band_mult = prev_band_mult,
  n_cores = 4
)

# After theta_ret <- polish_many(...)

# 1) build hazards per world (replace expand_hazards call with your actual helper)
haz_ret_list <- purrr::imap(theta_ret, ~ {
  expand_hazards(
    pars = .x, age = age, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    dplyr::select(age, trans, hazard) |>
    dplyr::mutate(world = .y, system = "returns")
})
source("R/simulation_functions2.R")
haz_snap_hu <- snap_hazards_agewise(
  haz_long_base  = haz_ret_list[[1]] |> dplyr::select(age, trans, hazard),
  ground_summary = ground_summary,
  age            = age,
  age_int        = age_int,
  trans          = c("hu","uh"),
  bounds         = c(1e-12, 2),
  smooth_w       = 0,     # set >0 to tame jitter
  verbose        = TRUE,
  maxit=1000
)

haz_snap_hu |> 
  mutate(version = "snapped",
         system="returns") |> 
  bind_rows(haz_ret |> filter(world==1)) |> 
  ggplot(aes(x=age,y=hazard,color=trans, linetype=version))+
  geom_line()+
  scale_y_log10()
lt_snap <- run_mslt_from_hazards(haz_snap_hu, age=age, init=init, age_int=age_int)

plot(lt_snap$prevalence - ground_summary$prevalence)
plot(lt_snap$lx - ground_summary$lx)








# 2) apply bump calibration to each returns world (HU/UH only)
cal_ret <- purrr::map(haz_ret_list, ~ calibrate_hazards_multibump(
  haz_long_base = dplyr::select(.x, age, trans, hazard),
  ground_summary = ground_summary,
  age = age,
  init = init,
  age_int = age_int,
  trans = c("hu","uh"),
  centers = seq(55, 90, by = 3),
  width = 4,
  w_int = 25,
  integrand = "health",
  w_prev = 0.5,
  w_lx = 0.1,
  ridge = 0.5,
  clamp_coef = c(-1, 1),
  clamp_haz = c(1e-12, 2),
  maxit = 200,
  param_mode = "balance",
  delta_only = F,     # <--- key
  ridge_gamma = 2
))
haz_ret_cal <- purrr::imap_dfr(cal_ret, ~ .x$haz_adj |> dplyr::mutate(world = .y, system="returns",version="post"))


# 3) now run mslt from hazards for diagnostics

lt_ret_cal <- purrr::imap_dfr(cal_ret, ~ {
  run_mslt_from_hazards(.x$haz_adj, age = age, init = init, age_int = age_int) |>
    dplyr::mutate(world = .y, system = "returns")
})
lt_ret_cal |> 
  select(world, age, lx) |> 
  pivot_wider(names_from = world, values_from = lx) |> 
  mutate(base_pattern = `5`) |> 
  pivot_longer(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`), names_to = "world", values_to = "lx") |> 
  mutate(residual = base_pattern - lx) |> 
  ggplot(aes(x = age,y=residual, color = world)) +
  geom_line() +
  labs(title = "lx residuals")

# check that hu or uh changed somewhere


hle_vec <- lt_ret_cal |>
  dplyr::group_by(world) |>
  summarize(hle = sum(lx * (1-prevalence))) |> 
  pull(hle)
hle_ground <- ground_summary |> summarize(hle = sum(lx*(1-prevalence))) |> pull(hle)
round(hle_vec - hle_ground,5)
haz_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  expand_hazards(
    th,
    age = age,
    age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(world = i, system = "returns",version="pre")
})
haz_ret |> bind_rows(haz_ret_cal) |> 
  filter(trans %in%c("hu","uh")) |> 
  pivot_wider(names_from = version,values_from = hazard) |> 
  mutate(ratio=post/pre) |> 
  ggplot(aes(x=age, y = ratio, color=trans))+
  geom_line()+
  scale_y_log10()+
  facet_wrap(~world)

compare_haz
lt_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  expand_hazards(
    th,
    age = age,
    age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(world = i, system = "returns")
}) |>
  group_split(world) |>
  map_dfr(function(hz) {
    w <- unique(hz$world)
    sys <- unique(hz$system)
    P <- haz_to_probs(hz |> select(age, trans, hazard), age = age, age_int = age_int)
    calculate_lt(P, init = init) |>
      mutate(world = w, system = sys)
  })
hle_pre <-
lt_ret |> 
  group_by(world) |> 
  summarize(hle = sum(lx * (1-prevalence))) |> 
  pull(hle)
round(hle_pre - hle_ground,5)
# ---------
# (3) RETURNS selection: keep good fits, then maximize diversity (returns only)
# ---------

fit_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  score_fit(
    th,
    ground_summary,
    age = age,
    init = init,
    age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(pert = i, system = "returns")
})

keep_ret <- fit_ret |>
  filter(rmse_lx <= rmse_lx_max, rmse_prev <= rmse_prev_max)

theta_pool <- theta_ret[keep_ret$pert]

# If we have too few, relax thresholds slightly (rare but can happen)
if (length(theta_pool) < (K_final - 1)) {
  warning("Too few returns worlds passed thresholds; using the best-fitting ones instead.", call. = FALSE)
  theta_pool <- theta_ret[order(fit_ret$rmse_prev + fit_ret$rmse_lx)][seq_len(min(length(theta_ret), K_final - 1))]
}

# Diversity selection among the pool (excluding the ground world, which we will force in)
H_pool <- do.call(rbind, lapply(theta_pool, haz_feat_from_theta))
D_pool <- as.matrix(dist(H_pool))

sel_pool <- select_farthest(D_pool, k = K_final - 1)
theta_sel <- theta_pool[sel_pool]

# Force the ground world into position world==5 for downstream tiny test
# (world numbering is purely for presentation / plotting)
theta_chosen <- vector("list", K_final)
theta_chosen[[5]] <- ground_pars_pw

# Fill remaining slots with selected returns worlds
fill_slots <- setdiff(seq_len(K_final), 5)
theta_chosen[fill_slots] <- theta_sel[seq_len(min(length(theta_sel), length(fill_slots)))]

# Hazards for chosen returns worlds (world id 1..K_final)
haz_chosen <- purrr::imap_dfr(theta_chosen, function(th, i) {
  expand_hazards(
    th,
    age = age,
    age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    mutate(world = i, system = "returns")
})

# ---------
# (4) NO-RETURNS worlds: plot Rx from returns, eyeball-pick two, derive hazards, swap
# ---------

Rx_plot_df <-
  haz_chosen |>
  filter(trans %in% c("hd","ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd) |>
  select(age, world, Rx)

ggplot(Rx_plot_df, aes(x = age, y = Rx, color = as.factor(world))) +
  geom_line() +
  labs(color = "returns world", y = "Rx = ud/hd",
       title = "Rx curves from chosen returns worlds (pick low/high)")


Rx_source_worlds <- c(3, 4)  # returns-world IDs to supply Rx
nr_world_ids     <- c(1, 6)  # which world IDs to replace with no-returns
stopifnot(length(Rx_source_worlds) == length(nr_world_ids))
stopifnot(length(Rx_source_worlds) == K_noreturn)

haz_nr_list <- purrr::map2(Rx_source_worlds, nr_world_ids, function(src_w, out_w) {
  
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

# Replace selected world IDs with their no-return versions
haz_chosen <-
  haz_chosen |>
  filter(!(world %in% nr_world_ids)) |>
  bind_rows(haz_nr) |>
  arrange(world, trans, age)

# (Optional) lifetables for no-return worlds (useful for checks)
haz_nr_list2 <- split(haz_nr |> select(age, trans, hazard), haz_nr$world)

lt_nr <- purrr::map(haz_nr_list2, function(hz) {
  P <- haz_to_probs(hz, age = age, age_int = age_int)
  calculate_lt(P, init = init)
})

# ---------
# (5) Export + quick plots
# ---------

dir.create("data", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(haz_chosen, "data/haz_octet.csv.gz")

haz_chosen |>
  filter(hazard > 1e-10) |> 
  ggplot(aes(x = age, y = hazard, color = as.factor(world))) +
  geom_line() +
  facet_wrap(~trans, scales = "free_y") +
  labs(color = "world", title = "Hazards for chosen worlds (returns + no-returns)") +
  scale_y_log10()

# ---------
# (6) Tiny sanity test: world 5 is the ground world
# ---------

lt_world5 <- calculate_lt(
  haz_to_probs(
    haz_chosen |> filter(world == 5) |> select(age, trans, hazard),
    age = age,
    age_int = age_int
  ),
  init = init
)

print(max(abs(lt_world5$lx - ground_summary$lx)))
print(max(abs(lt_world5$prevalence - ground_summary$prevalence)))

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
hle_sullivan <-lt_chosen |> 
  group_by(world) |> 
  summarize(hle_sullivan = sum((1-prevalence) * lx),
            hle_mslt = sum(lh)) |> 
  pull(hle_sullivan)
round(hle_sullivan - hle_sullivan[5],5)

lt_chosen |> 
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


