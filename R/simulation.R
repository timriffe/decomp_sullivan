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
# ---------
# (3) NO-RETURNS system: derive hazards directly from Rx picked by eyeballing
# ---------

# 1) Plot Rx(age) = ud/hd for the chosen RETURNS worlds (from Section 4)
Rx_plot_df <-
  haz_chosen |>
  dplyr::filter(system == "returns", trans %in% c("hd","ud")) |>
  tidyr::pivot_wider(names_from = trans, values_from = hazard) |>
  dplyr::mutate(Rx = ud / hd) |>
  dplyr::select(age, world, Rx)

ggplot2::ggplot(Rx_plot_df, ggplot2::aes(x = age, y = Rx, color = as.factor(world))) +
  ggplot2::geom_line() +
  ggplot2::labs(color = "returns world", y = "Rx = ud/hd")

# 2) YOU set these after eyeballing the plot:
#    - which returns worlds supply Rx (usually 2 worlds: low/high)
#    - which world IDs in the final set you want to become no-returns
Rx_source_worlds <- c(3, 4)  # <-- change after eyeballing
nr_world_ids     <- c(1, 6)  # <-- which worlds to replace with no-returns

stopifnot(length(Rx_source_worlds) == length(nr_world_ids))

# 3) Derive no-return hazards from those Rx curves (machine-precise match to lx & prev)
haz_nr_list <- purrr::map2(Rx_source_worlds, nr_world_ids, function(src_w, out_w) {
  
  Rx_vec <- Rx_plot_df |>
    dplyr::filter(world == src_w) |>
    dplyr::arrange(age) |>
    dplyr::pull(Rx)
  
  derive_noreturns_hazards(
    age  = ground_summary$age,
    lx   = ground_summary$lx,
    prev = ground_summary$prevalence,
    Rx   = Rx_vec,
    age_int = age_int,
    verbose = FALSE
  ) |>
    dplyr::mutate(world = out_w, system = "noreturn")
})

haz_nr <- dplyr::bind_rows(haz_nr_list)

# 4) Replace the selected world IDs with their no-return versions
#    (i.e., drop those returns worlds and bind no-returns in)
haz_chosen <-
  haz_chosen |>
  dplyr::filter(!(system == "returns" & world %in% nr_world_ids)) |>
  dplyr::bind_rows(haz_nr) |>
  dplyr::arrange(world, trans, age)

# optional quick check plot
haz_chosen |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = hazard, color = as.factor(world))) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~trans, scales = "free_y") +
  ggplot2::labs(color = "world")

# export
readr::write_csv(haz_chosen, "data/haz_octet.csv.gz")


# ---------
# (4) Combine, score fit, then select 5-8 diverse worlds

# Combine returns (thetas) + noreturns (hazards/lt)
n_ret <- length(theta_ret)
n_nr  <- length(haz_nr)

system_all <- c(rep("returns", n_ret), rep("noreturn", n_nr))

# score returns via existing score_fit()
fit_ret <- purrr::imap_dfr(theta_ret, function(th, i) {
  score_fit(th, ground_summary, age, init=init, age_int=age_int, knot_age=knot_age) |>
    mutate(pert = i, system = "returns")
})

# score noreturns directly from their lt objects
fit_nr <- purrr::imap_dfr(lt_nr, function(lt, j) {
  tibble::tibble(
    rmse_lx   = sqrt(mean((lt$lx - ground_summary$lx)^2)),
    rmse_prev = sqrt(mean((lt$prevalence - ground_summary$prevalence)^2)),
    maxabs_lx   = max(abs(lt$lx - ground_summary$lx)),
    maxabs_prev = max(abs(lt$prevalence - ground_summary$prevalence))
  ) |>
    mutate(pert = n_ret + j, system = "noreturn")
})

fit_tab <- dplyr::bind_rows(fit_ret, fit_nr)


# had to eyeball limits. No returns-systems are harder to match
keep <- fit_tab |>
  filter(rmse_lx <= rmse_lx_max, rmse_prev <= rmse_prev_max)


keep_ids    <- keep$pert
theta_keep  <- theta_all[keep_ids]
system_keep <- system_all[keep_ids]

# Diversity distance on stacked log hazards
# Build hazard feature matrix for both systems
haz_matrix_from_haz <- function(haz_df) {
  haz_df |>
    dplyr::mutate(loghaz = log(hazard)) |>
    dplyr::arrange(trans, age) |>
    dplyr::pull(loghaz)
}

# keep objects
keep_ids <- keep$pert

theta_keep <- theta_ret[keep_ids[keep_ids <= n_ret]]
nr_keep_ids <- keep_ids[keep_ids > n_ret] - n_ret

haz_keep <- c(
  lapply(theta_keep, function(th) expand_hazards(th, age=age, age_int=age_int, knot_age=knot_age)),
  haz_nr[nr_keep_ids]
)

system_keep <- c(
  rep("returns", length(theta_keep)),
  rep("noreturn", length(nr_keep_ids))
)

H <- do.call(rbind, lapply(haz_keep, haz_matrix_from_haz))
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
haz_chosen   <- haz_keep[sel]
sys_chosen    <- system_keep[sel]


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
  ggplot(aes(x=age,y=Ra,color = as.factor(world), linetype = system)) +
  geom_line(linewidth=1.2) +
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


# -------------------------------------- #
# redux for no-returns to get exact fits
# -------------------------------------- #
Ra_ssumptions <-
haz_chosen |> 
  filter(trans %in% c("hd","ud"),
         world %in% c(3,4)) |> 
  pivot_wider(names_from=trans,values_from=hazard) |> 
  mutate(Ra = ud / hd) 


haz1 <- derive_noreturns_hazards(ground_summary$age, 
                         lx = ground_summary$lx,
                         prev = ground_summary$prevalence,
                         Rx = Ra_ssumptions |> filter(world == 3) |> pull(Ra))|> 
  mutate(world = 1,
         system = "noreturn")


haz2 <- derive_noreturns_hazards(ground_summary$age, 
                                 lx = ground_summary$lx,
                                 prev = ground_summary$prevalence,
                                 Rx = Ra_ssumptions |> filter(world == 4) |> pull(Ra))|> 
  mutate(world = 6,
         system = "noreturn")

haz_chosen <-
  haz_chosen |> 
  filter(system == "returns") |> 
  bind_rows(haz1) |> 
  bind_rows(haz2)

write_csv(haz_chosen,"data/haz_octet.csv.gz")

haz_chosen |> 
  ggplot(aes(x=age,y=hazard,color = as.factor(world))) +
  geom_line() +
  facet_wrap(~trans)
