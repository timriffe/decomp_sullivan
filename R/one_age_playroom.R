# ------------------------------------------------------------
# One-age experimentation sandbox
# Robust to long-format hazard objects with extra grouping vars
# ------------------------------------------------------------

strip_num <- function(x) {
  as.numeric(unname(x))
}

hazards_to_sr <- function(hu, uh) {
  hu <- strip_num(hu)
  uh <- strip_num(uh)
  c(
    s = log(hu + uh),
    r = log(hu / uh)
  )
}

sr_to_hazards <- function(s, r) {
  c(
    hu = exp(s + r) / (1 + exp(r)),
    uh = exp(s)     / (1 + exp(r))
  )
}

state_to_lx_prev <- function(lh1, lu1) {
  lh1 <- strip_num(lh1)
  lu1 <- strip_num(lu1)
  lx1 <- lh1 + lu1
  prev <- if (is.na(lx1) || lx1 <= 0) NA_real_ else lu1 / lx1
  c(lx = lx1, prev = prev)
}

.one_age_step_local <- function(lh0, lu0, hd, ud, hu, uh, age_int) {
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("Package 'expm' is needed for the local one-age fallback.", call. = FALSE)
  }
  Q <- matrix(c(
    -(hu + hd),  hu,         hd,
    uh,        -(uh + ud),  ud,
    0,          0,          0
  ), nrow = 3, byrow = TRUE)
  
  P <- expm::expm(Q * age_int)
  start <- c(lh0, lu0, 1 - lh0 - lu0)
  out <- as.numeric(start %*% P)
  c(lh1 = out[1], lu1 = out[2], ld1 = out[3])
}

.step_one_age <- function(lh0, lu0, hd, ud, hu, uh, age_int) {
  if (exists("one_age_step", mode = "function")) {
    out <- one_age_step(
      lh0 = lh0, lu0 = lu0,
      hd = hd, ud = ud, hu = hu, uh = uh,
      age_int = age_int
    )
    return(c(
      lh1 = strip_num(out["lh1"]),
      lu1 = strip_num(out["lu1"])
    ))
  }
  .one_age_step_local(lh0, lu0, hd, ud, hu, uh, age_int)[c("lh1", "lu1")]
}

.extract_rx_vec <- function(world, age) {
  rx_df <- Rx_plot_df |>
    dplyr::filter(world == !!world) |>
    dplyr::select(age, Rx) |>
    dplyr::distinct() |>
    dplyr::arrange(age) |>
    dplyr::group_by(age) |>
    dplyr::summarise(Rx = dplyr::first(Rx), .groups = "drop")
  
  rx_match <- match(age, strip_num(rx_df$age))
  Rx_vec <- strip_num(rx_df$Rx[rx_match])
  
  if (length(Rx_vec) != length(age) || anyNA(Rx_vec)) {
    stop(
      "Could not build a full Rx vector aligned to ground_summary$age. ",
      "Check that Rx_plot_df has one Rx value per age for the requested world.",
      call. = FALSE
    )
  }
  
  Rx_vec
}

.extract_long_hazard_vec <- function(df, trans_name, age, world, prefer_system = NULL) {
  z <- df |>
    dplyr::filter(world == !!world, trans == !!trans_name)
  
  if (!is.null(prefer_system) && "system" %in% names(z)) {
    z1 <- z |> dplyr::filter(system == prefer_system)
    if (nrow(z1) > 0) z <- z1
  }
  
  z <- z |>
    dplyr::select(age, hazard) |>
    dplyr::distinct() |>
    dplyr::arrange(age) |>
    dplyr::group_by(age) |>
    dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
  
  idx <- match(age, strip_num(z$age))
  out <- strip_num(z$hazard[idx])
  
  if (length(out) != length(age) || anyNA(out)) {
    stop(
      "Could not build a full hazard vector for trans='", trans_name,
      "' aligned to ground_summary$age. Check haz_ret_pass1 structure.",
      call. = FALSE
    )
  }
  
  out
}

.extract_pass1_pair <- function(world, age) {
  if (!is.data.frame(haz_ret_pass1) ||
      !all(c("world", "trans", "age", "hazard") %in% names(haz_ret_pass1))) {
    stop("haz_ret_pass1 must be a long data frame with columns world, trans, age, hazard.",
         call. = FALSE)
  }
  
  list(
    hu0 = .extract_long_hazard_vec(haz_ret_pass1, "hu", age, world, prefer_system = NULL),
    uh0 = .extract_long_hazard_vec(haz_ret_pass1, "uh", age, world, prefer_system = NULL)
  )
}

.extract_snap_pair <- function(world, age, lx, prev, Rx_vec, age_int) {
  pass1 <- .extract_pass1_pair(world, age)
  
  snapped <- derive_returns_hazards_from_Rx(
    age = age,
    lx = lx,
    prev = prev,
    Rx = Rx_vec,
    age_int = age_int,
    hu_0 = pass1$hu0,
    uh_0 = pass1$uh0,
    verbose = FALSE
  )
  
  if (is.data.frame(snapped) && all(c("trans", "age", "hazard") %in% names(snapped))) {
    hu_snap <- snapped |>
      dplyr::filter(trans == "hu") |>
      dplyr::select(age, hazard) |>
      dplyr::distinct() |>
      dplyr::arrange(age) |>
      dplyr::group_by(age) |>
      dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
    uh_snap <- snapped |>
      dplyr::filter(trans == "uh") |>
      dplyr::select(age, hazard) |>
      dplyr::distinct() |>
      dplyr::arrange(age) |>
      dplyr::group_by(age) |>
      dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
    
    hu_snap <- strip_num(hu_snap$hazard[match(age, strip_num(hu_snap$age))])
    uh_snap <- strip_num(uh_snap$hazard[match(age, strip_num(uh_snap$age))])
  } else if (is.data.frame(snapped) && all(c("hu", "uh") %in% names(snapped))) {
    hu_snap <- strip_num(snapped$hu)
    uh_snap <- strip_num(snapped$uh)
  } else if (is.list(snapped) && all(c("hu", "uh") %in% names(snapped))) {
    hu_snap <- strip_num(snapped$hu)
    uh_snap <- strip_num(snapped$uh)
  } else {
    stop("Could not extract snapped hu/uh from derive_returns_hazards_from_Rx() output.",
         call. = FALSE)
  }
  
  if (length(hu_snap) != length(age) || length(uh_snap) != length(age) ||
      anyNA(hu_snap) || anyNA(uh_snap)) {
    stop("Snapped hu/uh could not be aligned to the full age vector.", call. = FALSE)
  }
  
  list(
    hu0 = pass1$hu0,
    uh0 = pass1$uh0,
    hu_snap = hu_snap,
    uh_snap = uh_snap
  )
}

eval_sr_one_age <- function(play, s, r) {
  hz <- sr_to_hazards(s, r)
  hu <- strip_num(hz["hu"])
  uh <- strip_num(hz["uh"])
  
  pred <- .step_one_age(
    lh0 = play$lh0,
    lu0 = play$lu0,
    hd = play$hd,
    ud = play$ud,
    hu = hu,
    uh = uh,
    age_int = play$age_int
  )
  
  lh1 <- strip_num(pred["lh1"])
  lu1 <- strip_num(pred["lu1"])
  st <- state_to_lx_prev(lh1, lu1)
  
  data.frame(
    s = strip_num(s),
    r = strip_num(r),
    hu = hu,
    uh = uh,
    lh1 = lh1,
    lu1 = lu1,
    lx1 = strip_num(st["lx"]),
    prev1 = strip_num(st["prev"]),
    lh_resid = lh1 - play$lh_target,
    lu_resid = lu1 - play$lu_target,
    lx_resid = strip_num(st["lx"]) - play$lx_target,
    prev_resid = strip_num(st["prev"]) - play$prev_target,
    sq_loss_state = (lh1 - play$lh_target)^2 + (lu1 - play$lu_target)^2
  )
}

optimize_r_given_s <- function(play, s, r_window = log(3), tol = 1e-14) {
  r0 <- play$r_snap
  f <- function(r) eval_sr_one_age(play, s, r)$sq_loss_state
  
  opt <- optimize(
    f,
    interval = c(r0 - r_window, r0 + r_window),
    tol = tol
  )
  
  eval_sr_one_age(play, s, opt$minimum)
}

sweep_s_optimize_r <- function(play, s_values, r_window = log(3), tol = 1e-14) {
  out <- lapply(s_values, function(s) optimize_r_given_s(play, s, r_window, tol))
  do.call(rbind, out)
}

plot_sweep_s <- function(sw, y = "sq_loss_state") {
  plot(sw$s, sw[[y]], type = "l", xlab = "s", ylab = y)
}

optim_sr_local <- function(play, objective = c("state", "lxprev"),
                           s_window = log(1.5), r_window = log(3),
                           maxit = 1000) {
  objective <- match.arg(objective)
  
  par0 <- c(play$s_snap, play$r_snap)
  lower <- c(play$s_snap - s_window, play$r_snap - r_window)
  upper <- c(play$s_snap + s_window, play$r_snap + r_window)
  
  fn <- function(par) {
    ev <- eval_sr_one_age(play, par[1], par[2])
    if (objective == "state") {
      ev$sq_loss_state
    } else {
      ev$lx_resid^2 + ev$prev_resid^2
    }
  }
  
  opt <- optim(
    par = par0,
    fn = fn,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = maxit, factr = 1, pgtol = 0)
  )
  
  list(
    opt = opt,
    eval = eval_sr_one_age(play, opt$par[1], opt$par[2])
  )
}

setup_one_age_playroom <- function(age_pick = 70, world = 1, verbose = TRUE) {
  age <- strip_num(ground_summary$age)
  lx <- strip_num(ground_summary$lx)
  prev <- strip_num(ground_summary$prevalence)
  
  idx <- which(age == age_pick)
  if (length(idx) != 1L) {
    stop("age_pick must match exactly one age in ground_summary$age.", call. = FALSE)
  }
  if (idx >= length(age)) {
    stop("age_pick must not be the terminal/open age.", call. = FALSE)
  }
  
  Rx_vec <- .extract_rx_vec(world = world, age = age)
  
  nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev = prev,
    Rx = Rx_vec,
    age_int = age_int,
    verbose = FALSE
  )
  
  if (is.data.frame(nr) && all(c("trans", "age", "hazard") %in% names(nr))) {
    hd <- nr |> dplyr::filter(trans == "hd") |> dplyr::select(age, hazard) |> dplyr::distinct() |> dplyr::arrange(age) |> dplyr::group_by(age) |> dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
    ud <- nr |> dplyr::filter(trans == "ud") |> dplyr::select(age, hazard) |> dplyr::distinct() |> dplyr::arrange(age) |> dplyr::group_by(age) |> dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
    hu_nr <- nr |> dplyr::filter(trans == "hu") |> dplyr::select(age, hazard) |> dplyr::distinct() |> dplyr::arrange(age) |> dplyr::group_by(age) |> dplyr::summarise(hazard = dplyr::first(hazard), .groups = "drop")
    
    hd <- strip_num(hd$hazard[match(age, strip_num(hd$age))])
    ud <- strip_num(ud$hazard[match(age, strip_num(ud$age))])
    hu_nr <- strip_num(hu_nr$hazard[match(age, strip_num(hu_nr$age))])
  } else if (is.data.frame(nr) && all(c("hd", "ud", "hu") %in% names(nr))) {
    hd <- strip_num(nr$hd)
    ud <- strip_num(nr$ud)
    hu_nr <- strip_num(nr$hu)
  } else if (is.list(nr) && all(c("hd", "ud", "hu") %in% names(nr))) {
    hd <- strip_num(nr$hd)
    ud <- strip_num(nr$ud)
    hu_nr <- strip_num(nr$hu)
  } else {
    stop("Could not extract hd/ud/hu from derive_noreturns_hazards() output.", call. = FALSE)
  }
  
  snap <- .extract_snap_pair(world, age, lx, prev, Rx_vec, age_int)
  
  hu_snap_i <- snap$hu_snap[idx]
  uh_snap_i <- snap$uh_snap[idx]
  sr <- hazards_to_sr(hu_snap_i, uh_snap_i)
  
  lh0 <- lx[idx] * (1 - prev[idx])
  lu0 <- lx[idx] * prev[idx]
  lh_target <- lx[idx + 1] * (1 - prev[idx + 1])
  lu_target <- lx[idx + 1] * prev[idx + 1]
  
  play <- list(
    age = age_pick,
    world = world,
    idx = idx,
    age_int = age_int,
    Rx_vec = Rx_vec,
    hd = hd[idx],
    ud = ud[idx],
    hu_nr = hu_nr[idx],
    hu0 = snap$hu0[idx],
    uh0 = snap$uh0[idx],
    hu_snap = hu_snap_i,
    uh_snap = uh_snap_i,
    s_snap = strip_num(sr["s"]),
    r_snap = strip_num(sr["r"]),
    lh0 = lh0,
    lu0 = lu0,
    lh_target = lh_target,
    lu_target = lu_target,
    lx_target = lx[idx + 1],
    prev_target = prev[idx + 1]
  )
  
  if (verbose) {
    test <- eval_sr_one_age(play, play$s_snap, play$r_snap)
    cat("one-age playroom set for age", age_pick, ", world", world, "\n", sep = " ")
    cat("  hd =", format(play$hd, digits = 12),
        ", ud =", format(play$ud, digits = 12), "\n", sep = "")
    cat("  no-returns hu =", format(play$hu_nr, digits = 12), "\n", sep = "")
    cat("  pass1 guess    hu =", format(play$hu0, digits = 12),
        ", uh =", format(play$uh0, digits = 12), "\n", sep = "")
    cat("  snapped pair   hu =", format(play$hu_snap, digits = 12),
        ", uh =", format(play$uh_snap, digits = 12), "\n", sep = "")
    cat("  s =", format(play$s_snap, digits = 12),
        ", r =", format(play$r_snap, digits = 12), "\n", sep = "")
    cat("  loss at snapped pair =", format(test$sq_loss_state, digits = 12), "\n", sep = "")
    cat("  lx residual at snapped pair   =",
        format(test$lx_resid, digits = 12), "\n", sep = "")
    cat("  prev residual at snapped pair =",
        format(test$prev_resid, digits = 12), "\n", sep = "")
  }
  
  play
}
play70 <- setup_one_age_playroom(age_pick = 70, world = 1, verbose = TRUE)
play75 <- setup_one_age_playroom(age_pick = 75, world = 1, verbose = TRUE)
play80 <- setup_one_age_playroom(age_pick = 80, world = 1, verbose = TRUE)
sw70 <- sweep_s_optimize_r(
  play70,
  s_values = seq(log(play70$hu_nr),log(play70$hu_nr)+5,length=1000),
  r_window = 100,
  tol = 1e-20
)
sw75 <- sweep_s_optimize_r(
  play75,
  s_values = seq(log(play75$hu_nr),log(play75$hu_nr)+5,length=1000),
  r_window = 100,
  tol = 1e-20
)
sw80 <- sweep_s_optimize_r(
  play80,
  s_values = seq(log(play80$hu_nr),log(play80$hu_nr)+5,length=1000),
  r_window = 100,
  tol = 1e-20
)

bind_rows(sw70 |> mutate(age=70),
          sw75 |> mutate(age=75),
          sw80 |> mutate(age=80)) |> 
  mutate(turnover = uh+hu) |> 
  ggplot(aes(x =turnover,y=prev_resid,color=as.factor(age)))+
  geom_line()
