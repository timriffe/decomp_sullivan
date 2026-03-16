# Drop-in replacement for derive_returns_hazards_from_Rx()
# Fresh-eye version following the intended 5-step routine:
#  1) derive hd, ud, hu_nr from no-returns world
#  2) fix p1 = (hu_nr, 0)
#  3) define p2 by perturbing hu upward and solving uh
#  4) orthogonally project (hu0, uh0) onto the detected line
#  5) with hu fixed at the projected value, optimize uh one last time
#
# Source AFTER simulation_functions2.R

derive_returns_hazards_from_Rx <- function(
    age,
    lx,
    prev,
    Rx,
    age_int = 1,
    bounds = c(1e-12, 2),
    hu_bounds = NULL,
    uh_bounds = NULL,
    hu_0,
    uh_0,
    turnover_K = 1.2,
    snap_tol = 1e-20,
    verbose = FALSE,
    eps_log = 1e-2,
    tiny = 1e-14,
    line_tol = 1e-12,
    final_uh_refine = TRUE,
    final_uh_tol = 1e-20,
    final_uh_accept_only_if_improves = TRUE,
    final_uh_skip_last_age = TRUE,
    extrap_last_age = TRUE,
    prevalence_point = NULL
){
  age  <- unname(as.numeric(age))
  lx   <- unname(as.numeric(lx))
  prev <- unname(as.numeric(prev))
  Rx   <- unname(as.numeric(Rx))
  hu_0 <- unname(as.numeric(hu_0))
  uh_0 <- unname(as.numeric(uh_0))
  
  n <- length(age)
  stopifnot(
    length(lx) == n,
    length(prev) == n,
    length(Rx) == n,
    length(hu_0) == n,
    length(uh_0) == n
  )
  
  # 1) fixed no-returns backbone: hd, ud, hu_nr
  hz_nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev = prev,   # intentionally passed through as supplied
    Rx = Rx,
    age_int = age_int,
    verbose = FALSE
  )
  
  if (!all(c("age", "trans", "hazard") %in% names(hz_nr))) {
    stop("derive_noreturns_hazards() must return columns age, trans, hazard.", call. = FALSE)
  }
  
  hd    <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
  ud    <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
  hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
  
  if (length(hd) != n || length(ud) != n || length(hu_nr) != n) {
    stop("No-returns hazards could not be aligned by age.", call. = FALSE)
  }
  
  if (is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  if (is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  
  hu_lo <- pmax(hu_bounds[, 1], bounds[1])
  hu_hi <- hu_bounds[, 2]
  uh_lo <- pmax(0, uh_bounds[, 1])
  uh_hi <- uh_bounds[, 2]
  
  # point-prevalence state targets for exact-age lh/lu reconstruction
  states <- build_states_targets(
    age = age,
    lx = lx,
    prev = prev,
    prevalence_point = prevalence_point
  )
  
  hu_out <- numeric(n)
  uh_out <- numeric(n)
  
  p2_loss <- rep(NA_real_, n)
  proj_loss <- rep(NA_real_, n)
  post_final_loss <- rep(NA_real_, n)
  accepted_final_uh <- rep(FALSE, n)
  degenerate <- rep(FALSE, n)
  
  for (i in seq_len(n)) {
    cap <- turnover_K * (hu_0[i] + uh_0[i])
    
    # 2) fixed anchor from no-returns world
    p1 <- c(clamp(hu_nr[i], hu_lo[i], hu_hi[i]), 0)
    p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
    
    # incoming pass-1 point
    p0 <- c(
      clamp(hu_0[i], hu_lo[i], hu_hi[i]),
      clamp(uh_0[i], uh_lo[i], uh_hi[i])
    )
    
    # 3) second anchor by nudging hu upward and solving uh
    hu2_try <- clamp(p1[1] * exp(eps_log), hu_lo[i], hu_hi[i])
    uh_hi_eff <- min(uh_hi[i], cap - hu2_try)
    if (!is.finite(uh_hi_eff)) uh_hi_eff <- uh_hi[i]
    uh_hi_eff <- max(uh_hi_eff, uh_lo[i])
    
    sol2 <- solve_uh_given_hu(
      i, states, hd, ud,
      hu = hu2_try,
      uh_lo = uh_lo[i],
      uh_hi = uh_hi_eff,
      age_int = age_int,
      tiny = tiny
    )
    
    have_p2 <- is.finite(sol2$loss) &&
      is.finite(sol2$uh) &&
      (abs(hu2_try - p1[1]) > 0) &&
      (sol2$loss <= line_tol)
    
    if (!have_p2) {
      # fallback: nudge uh upward and solve hu
      uh_seed <- clamp(max(uh_lo[i], tiny) * exp(eps_log), uh_lo[i], uh_hi[i])
      hu_hi_eff <- min(hu_hi[i], cap - uh_seed)
      if (!is.finite(hu_hi_eff)) hu_hi_eff <- hu_hi[i]
      hu_hi_eff <- max(hu_hi_eff, hu_lo[i])
      
      sol2b <- solve_hu_given_uh(
        i, states, hd, ud,
        uh = uh_seed,
        hu_lo = hu_lo[i],
        hu_hi = hu_hi_eff,
        age_int = age_int
      )
      
      have_p2 <- is.finite(sol2b$loss) &&
        is.finite(sol2b$hu) &&
        (abs(uh_seed - p1[2]) > 0) &&
        (sol2b$loss <= line_tol)
      
      if (have_p2) {
        p2 <- c(sol2b$hu, uh_seed)
        p2_loss[i] <- sol2b$loss
      } else {
        # last-resort use whichever candidate is finite
        if (is.finite(sol2$loss) && is.finite(sol2$uh) && abs(hu2_try - p1[1]) > 0) {
          p2 <- c(hu2_try, sol2$uh)
          p2_loss[i] <- sol2$loss
          have_p2 <- TRUE
        } else if (is.finite(sol2b$loss) && is.finite(sol2b$hu) && abs(uh_seed - p1[2]) > 0) {
          p2 <- c(sol2b$hu, uh_seed)
          p2_loss[i] <- sol2b$loss
          have_p2 <- TRUE
        } else {
          p2 <- p1
          p2_loss[i] <- NA_real_
          have_p2 <- FALSE
        }
      }
    } else {
      p2 <- c(hu2_try, sol2$uh)
      p2_loss[i] <- sol2$loss
    }
    
    v <- p2 - p1
    
    if (sum(v * v) < 1e-30) {
      degenerate[i] <- TRUE
      hu_proj <- p1[1]
      uh_proj <- p1[2]
    } else {
      feas <- clamp_t_feasible(
        a = p1, v = v,
        hu_lo = hu_lo[i], hu_hi = hu_hi[i],
        uh_lo = uh_lo[i], uh_hi = uh_hi[i],
        cap = cap
      )
      
      if (!isTRUE(feas$ok)) {
        degenerate[i] <- TRUE
        hu_proj <- p1[1]
        uh_proj <- p1[2]
      } else {
        # 4) pure orthogonal projection onto the feasible line segment
        pr <- project_point_to_line(p0, p1, v)
        t_proj <- clamp(pr$t, feas$t_lo, feas$t_hi)
        p_proj <- p1 + t_proj * v
        hu_proj <- p_proj[1]
        uh_proj <- p_proj[2]
      }
    }
    
    proj_loss[i] <- loss_one_age(
      i, states, hd, ud,
      hu = hu_proj, uh = uh_proj,
      age_int = age_int
    )
    
    hu_hat <- hu_proj
    uh_hat <- uh_proj
    post_final_loss[i] <- proj_loss[i]
    
    # 5) one final uh optimization with hu fixed at projected hu
    do_final <- isTRUE(final_uh_refine) &&
      !(isTRUE(final_uh_skip_last_age) && i == n)
    
    if (do_final) {
      uh_hi_eff2 <- min(uh_hi[i], cap - hu_hat)
      if (!is.finite(uh_hi_eff2)) uh_hi_eff2 <- uh_hi[i]
      uh_hi_eff2 <- max(uh_hi_eff2, uh_lo[i])
      
      sol_final <- solve_uh_given_hu(
        i, states, hd, ud,
        hu = hu_hat,
        uh_lo = uh_lo[i],
        uh_hi = uh_hi_eff2,
        age_int = age_int,
        tiny = tiny
      )
      
      if (is.finite(sol_final$loss) && is.finite(sol_final$uh)) {
        improve <- proj_loss[i] - sol_final$loss
        accept <- if (isTRUE(final_uh_accept_only_if_improves)) {
          improve > 0
        } else {
          TRUE
        }
        
        if (accept) {
          uh_hat <- clamp(sol_final$uh, uh_lo[i], uh_hi_eff2)
          post_final_loss[i] <- sol_final$loss
          accepted_final_uh[i] <- TRUE
        }
      }
    }
    
    hu_out[i] <- hu_hat
    uh_out[i] <- uh_hat
    
    if (verbose && (i %% 10 == 0 || i == 1 || i == n)) {
      message(sprintf(
        "age %s: proj_loss=%.3e post_final=%.3e p2_loss=%s acc=%s deg=%s",
        age[i],
        proj_loss[i],
        post_final_loss[i],
        format(p2_loss[i], scientific = TRUE, digits = 2),
        accepted_final_uh[i],
        degenerate[i]
      ))
    }
  }
  
  if (isTRUE(extrap_last_age) && n >= 3) {
    min_haz <- bounds[1]
    extrap_last <- function(x) {
      if (is.finite(x[n - 1]) && is.finite(x[n - 2]) &&
          x[n - 1] > 0 && x[n - 2] > 0) {
        exp(2 * log(x[n - 1]) - log(x[n - 2]))
      } else if (is.finite(x[n - 1])) {
        x[n - 1]
      } else {
        min_haz
      }
    }
    
    hu_out[n] <- clamp(pmax(extrap_last(hu_out), min_haz), hu_lo[n], hu_hi[n])
    uh_out[n] <- clamp(pmax(extrap_last(uh_out), min_haz), uh_lo[n], uh_hi[n])
    
    cap_n <- turnover_K * (hu_0[n] + uh_0[n])
    if (is.finite(cap_n) && (hu_out[n] + uh_out[n] > cap_n)) {
      uh_out[n] <- max(uh_lo[n], cap_n - hu_out[n])
    }
  }
  
  out <- tibble::tibble(
    age = rep(age, times = 4),
    trans = rep(c("hd", "hu", "ud", "uh"), each = n),
    hazard = c(hd, hu_out, ud, uh_out)
  )
  
  attr(out, "diagnostics") <- list(
    p2_loss = p2_loss,
    proj_loss = proj_loss,
    post_final_loss = post_final_loss,
    accepted_final_uh = accepted_final_uh,
    degenerate = degenerate
  )
  
  out
}