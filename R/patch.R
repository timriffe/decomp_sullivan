# Drop-in replacement for derive_returns_hazards_from_Rx()
# Patch: recompute the no-returns anchor p1 on the SAME point-state target scale
# used by the returns loss, instead of importing hu_nr directly from the
# interval-prevalence no-returns solve.
#
# Source AFTER simulation_functions2.R

.recompute_p1_hu <- function(i, states, hd, ud, hu_lo, hu_hi, age_int,
                             tol = 1e-14, tiny = 1e-14) {
  f <- function(hu) {
    loss_one_age(
      i = i, states = states, hd = hd, ud = ud,
      hu = hu, uh = 0, age_int = age_int
    )
  }
  
  opt <- optimize(f, interval = c(hu_lo, hu_hi), tol = tol)
  hu_hat <- as.numeric(opt$minimum)
  ls_hat <- as.numeric(opt$objective)
  
  list(
    hu = hu_hat,
    loss = ls_hat
  )
}

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
    line_tol = 1e-12,
    eps_log = 1e-2,
    tiny = 1e-14,
    refine_tol = 1e-16,
    refine_dist_w = 1e-6,
    extrap_last_age = TRUE,
    verbose = FALSE,
    prevalence_point = NULL
) {
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
  
  hz_nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev = prev,
    Rx = Rx,
    age_int = age_int,
    verbose = FALSE
  )
  
  if (!all(c("age", "trans", "hazard") %in% names(hz_nr))) {
    stop("derive_noreturns_hazards() must return columns age, trans, hazard.",
         call. = FALSE)
  }
  
  hd <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
  ud <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
  hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
  
  if (is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  if (is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  
  hu_lo <- pmax(hu_bounds[, 1], bounds[1])
  hu_hi <- hu_bounds[, 2]
  uh_lo <- pmax(0, uh_bounds[, 1])
  uh_hi <- uh_bounds[, 2]
  
  states <- build_states_targets(
    age = age,
    lx = lx,
    prev = prev,
    prevalence_point = prevalence_point
  )
  
  hu_out <- numeric(n)
  uh_out <- numeric(n)
  
  p1_loss <- rep(NA_real_, n)
  p2_loss <- rep(NA_real_, n)
  degenerate <- rep(FALSE, n)
  
  for (i in seq_len(n)) {
    cap <- turnover_K * (hu_0[i] + uh_0[i])
    
    # PATCH: p1 is recomputed from the point-state target loss with uh = 0,
    # using the fixed hd, ud inherited from the interval-prevalence no-returns solve.
    p1_sol <- .recompute_p1_hu(
      i = i,
      states = states,
      hd = hd,
      ud = ud,
      hu_lo = hu_lo[i],
      hu_hi = hu_hi[i],
      age_int = age_int,
      tol = min(refine_tol, 1e-14),
      tiny = tiny
    )
    
    p1 <- c(clamp(p1_sol$hu, hu_lo[i], hu_hi[i]), 0)
    p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
    p1_loss[i] <- p1_sol$loss
    
    p0 <- c(
      clamp(hu_0[i], hu_lo[i], hu_hi[i]),
      clamp(uh_0[i], uh_lo[i], uh_hi[i])
    )
    
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
      hu_hat <- p1[1]
      uh_hat <- p1[2]
    } else {
      feas <- clamp_t_feasible(
        a = p1, v = v,
        hu_lo = hu_lo[i], hu_hi = hu_hi[i],
        uh_lo = uh_lo[i], uh_hi = uh_hi[i],
        cap = cap
      )
      
      if (!isTRUE(feas$ok)) {
        hu_hat <- p1[1]
        uh_hat <- p1[2]
      } else {
        pr <- project_point_to_line(p0, p1, v)
        t_proj <- clamp(pr$t, feas$t_lo, feas$t_hi)
        
        f_t <- function(t) {
          pt <- p1 + t * v
          ls <- loss_one_age(
            i, states, hd, ud,
            hu = pt[1], uh = pt[2],
            age_int = age_int
          )
          ls + refine_dist_w * (t - t_proj)^2
        }
        
        opt_t <- optimize(f_t, interval = c(feas$t_lo, feas$t_hi), tol = refine_tol)
        t_star <- opt_t$minimum
        
        p_star <- p1 + t_star * v
        hu_hat <- p_star[1]
        uh_hat <- p_star[2]
      }
    }
    
    hu_out[i] <- hu_hat
    uh_out[i] <- uh_hat
    
    if (verbose && (i %% 10 == 0 || i == 1 || i == n)) {
      ls <- loss_one_age(i, states, hd, ud, hu = hu_out[i], uh = uh_out[i], age_int = age_int)
      message(sprintf(
        "age %s: loss=%.3e  p1_loss=%s  p2_loss=%s  deg=%s",
        age[i],
        ls,
        format(p1_loss[i], scientific = TRUE, digits = 2),
        format(p2_loss[i], scientific = TRUE, digits = 2),
        degenerate[i]
      ))
    }
  }
  
  if (isTRUE(extrap_last_age) && n >= 3) {
    min_haz <- bounds[1]
    extrap_last <- function(x) {
      if (is.finite(x[n - 1]) && is.finite(x[n - 2]) && x[n - 1] > 0 && x[n - 2] > 0) {
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
    p1_loss = p1_loss,
    p2_loss = p2_loss,
    degenerate = degenerate,
    hu_nr = hu_nr
  )
  
  out
}