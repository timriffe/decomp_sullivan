
source("R/simulation_functions.R")

# -----------------------------------------------------------------------------
# top-level settings
# -----------------------------------------------------------------------------
age_int <- 1 #1/12
age     <- seq(50,100,by=age_int)
init    <- c(H = 1, U = 0)

# Hazard family for ground + candidates + polishing
haz_shape   <- "softkink"   # "loglinear", "piecewise", "softkink"
pivot_age   <- 70
shape_width <- 6

# Candidate stepping along null directions 
step_grid <- c(-0.3, -.1, .1, 0.3)   


# ---------------------------------------------------
# Ground world
# ---------------------------------------------------
# This arbitrary parameterization creates the ground world.
# from here onward we try to match the implied sullivan stats
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

ground_summary <- ground_mslt |> select(age, lx, prevalence, prevalence_point, Lx)

# Quick plot of targets
ground_summary |> select(-prevalence_point) |> 
  pivot_longer(c(lx, prevalence), names_to = "measure", values_to = "value") |>
  ggplot(aes(age, value, color = measure)) +
  geom_line() +
  labs(title = "Ground Sullivan targets")

# ---------------------------------------
# Null directions from Jacobian, creates world diversity
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
# ------------------------------------------------------------
# Candidate hazards directly from cand
# ------------------------------------------------------------
cand <- make_candidates(
  base_pars  = ground_pars,
  free_names = free_names,
  V2 = V2,
  grid_steps = step_grid,
  exclude_zero = TRUE
)

# candidate parameter vectors
theta_list <- cand$theta0

haz_cand <- purrr::imap_dfr(theta_list, function(th, i) {
  pars_i <- if (length(th) == 12) th else as_piecewise_pars(th)
  
  expand_hazards(
    pars_i,
    age = age,
    age_int = age_int,
    haz_shape = haz_shape,
    pivot_age = pivot_age,
    shape_width = shape_width
  ) |>
    mutate(world = i, system = "cand_returns")
})

# ------------------------------------------------------------
# Rx curves from candidates: these determine the worlds!
# ------------------------------------------------------------
Rx_cand_df <- haz_cand |>
  filter(trans %in% c("hd", "ud")) |>
  pivot_wider(names_from = trans, values_from = hazard) |>
  mutate(Rx = ud / hd) |>
  select(world, age, Rx)

ggplot(Rx_cand_df, aes(x = age, y = Rx, color = factor(world))) +
  geom_line() +
  labs(title = "Rx curves from null-direction candidates", color = "world")

# --------------------------------
# derive returns worlds, can take long time,
# here we use Rx, we do fancy dancing to get
# starting values, and then we do more fancy dancing
# to snap to precision
# --------------------------------

n_cores <- parallel::detectCores()

rx_groups <- Rx_cand_df %>%
  group_by(world) %>%
  group_split()

rx_keys <- Rx_cand_df %>%
  group_by(world) %>%
  group_keys()

old_plan <- future::plan()

haz_returns_exact <- tryCatch({
  
  future::plan(future::multisession, workers = n_cores)
  
  furrr::future_map2_dfr(
    .x = rx_groups,
    .y = rx_keys$world,
    .f = function(.x, .world) {
      
      Rx_vec <- .x %>%
        arrange(age) %>%
        pull(Rx)
      
      derive_returns_hazards(
        age = ground_summary$age,
        lx = ground_summary$lx,
        prev_pt = ground_summary$prevalence_point,
        Rx = Rx_vec,
        theta_start = theta_list[[.world]],
        age_int = age_int,
        pivot_age = 75,
        shape_width = 5,
        anchor_age = 70,
        correction_max_age = 55,
        min_haz = 1e-12,
        tol = 1e-20,
        maxiter = 2000,
        factr = 1e6,
        w_lx = 1,
        w_prev = 1,
        slope_floor = 1e-6,
        hd_bracket = c(1e-12, 5),
        uh_bracket = c(1e-12, 5),
        diag_tol = 1e-14,
        hu_bounds_mult = c(0.7, 3.0),
        n_bracket = 21,
        n_refine = 20,
        eps = 0.01,
        world_id = .world,
        verbose = TRUE
      ) %>%
        mutate(
          world = .world,
          system = "returns"
        )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  
}, finally = {
  future::plan(old_plan)
})

# ------------------------------------
# derive noreturns worlds, also based on just Rx and the ground_summary
# ------------------------------------
haz_noreturns <- Rx_cand_df %>%
  group_by(.data$world) %>%
  group_modify(~{
    w <- unique(.y$world)
    Rx_vec <- .x %>% arrange(age) %>% pull(Rx)
    
    derive_noreturns_hazards(
      age  = ground_summary$age,
      lx   = ground_summary$lx,
      prev_pt = ground_summary$prevalence_point,
      Rx   = Rx_vec,
      age_int = age_int,
      tol = 1e-20, 
      maxiter = 100, 
      hu_bracket = c(1e-12, 5),
      hd_bracket = c(1e-12, 5),
      verbose = FALSE
    ) %>%
      mutate(system = "noreturns")
  }) %>%
  ungroup()


# ----------------------
# Combine and calc lifetables
# ---------------------
haz_all <-
 bind_rows(haz_returns_exact, haz_noreturns)

lt_all <-
  haz_all %>%
  group_by(system, world) %>%
  group_modify(~{
    haz_block_to_lt(.x, init = init, age_int = age_int)
  }) %>%
  ungroup()

# check precision
hle <-
lt_all |> 
  group_by(system,world) |> 
  summarise(hle=sum(Lx*(1-prevalence_interval)),.groups="drop") |> 
  pull(hle)

print(hle,digits=12)
# ------------------------------------------------------------
# save outputs
# ------------------------------------------------------------
save_systems <- c("returns","noreturns")
if (age_int == 1){
  readr::write_csv(haz_all |> filter(system %in% save_systems), file = "data/worlds_hazards_annual.csv.gz")
  readr::write_csv(lt_all |> filter(system %in% save_systems), file = "data/worlds_mslt_annual.csv.gz")
}
if (age_int == 1/12){
  readr::write_csv(haz_all |> filter(system %in% save_systems), file = "data/worlds_hazards_monthly.csv.gz")
  readr::write_csv(lt_all |> filter(system %in% save_systems), file = "data/worlds_mslt_monthly.csv.gz")
}


