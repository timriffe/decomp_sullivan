# Experimental override for derive_noreturns_hazards()
# Uses interval-average unhealthy share (trapezoid proxy) instead of point prevalence
# for constructing the target unhealthy stock vector inside the no-returns solver.
#
# Source AFTER simulation_functions.R

interval_prev_trapezoid <- function(lx, prev, last = c("point", "carry")) {
  last <- match.arg(last)
  lx <- as.numeric(lx)
  prev <- as.numeric(prev)
  
  if (length(lx) != length(prev)) {
    stop("lx and prev must have the same length.", call. = FALSE)
  }
  n <- length(lx)
  if (n < 2) return(prev)
  
  lu <- lx * prev
  
  out <- numeric(n)
  out[1:(n - 1)] <- (lu[1:(n - 1)] + lu[2:n]) / (lx[1:(n - 1)] + lx[2:n])
  
  out[n] <- if (last == "point") prev[n] else out[n - 1]
  out
}

derive_noreturns_hazards <- function(age,
                                     lx,
                                     prev,
                                     Rx,
                                     age_int = 1,
                                     min_haz = 1e-12,
                                     # root finding controls
                                     tol = 1e-20,
                                     maxiter = 100,
                                     # brackets (expanded automatically if needed)
                                     hu_bracket = c(1e-12, 5),
                                     hd_bracket = c(1e-12, 5),
                                     verbose = TRUE,
                                     prev_mode = c("interval", "point"),
                                     interval_last = c("carry", "point")) {
  
  prev_mode <- match.arg(prev_mode)
  interval_last <- match.arg(interval_last)
  
  age <- as.numeric(age)
  n <- length(age)
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  
  lx <- as.numeric(lx)
  if (length(lx) != n) stop("lx must have same length as age.", call. = FALSE)
  
  prev <- as.numeric(prev)
  if (length(prev) != n) stop("prev must have same length as age.", call. = FALSE)
  
  if (length(Rx) == 1) Rx <- rep(as.numeric(Rx), n)
  Rx <- as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or length(age).", call. = FALSE)
  
  prev_used <- if (prev_mode == "interval") {
    interval_prev_trapezoid(lx = lx, prev = prev, last = interval_last)
  } else {
    prev
  }
  
  # target state occupancies
  lu_tgt <- prev_used * lx
  lh_tgt <- lx - lu_tgt
  
  # numeric safety
  Rx <- pmax(Rx, min_haz)
  lx <- pmax(lx, min_haz)
  lu_tgt <- pmin(pmax(lu_tgt, 0), lx)
  lh_tgt <- pmax(lh_tgt, 0)
  
  make_P <- function(hu, hd, Rx_i, dt) {
    ud <- Rx_i * hd
    uh <- min_haz
    Q <- matrix(c(
      -(hu + hd),  hu,  hd,
      uh, -(uh + ud), ud,
      0, 0, 0
    ), nrow = 3, byrow = TRUE)
    expm::expm(Q * dt)
  }
  
  step_forward <- function(H, U, P) {
    v_next <- c(H, U, 0) %*% P
    c(H = as.numeric(v_next[1]), U = as.numeric(v_next[2]))
  }
  
  hu <- rep(NA_real_, n)
  hd <- rep(NA_real_, n)
  ud <- rep(NA_real_, n)
  uh <- rep(min_haz, n)
  
  for (i in 1:(n - 1)) {
    H0 <- lh_tgt[i]
    U0 <- lu_tgt[i]
    
    if ((H0 + U0) <= min_haz) {
      hu[i] <- min_haz
      hd[i] <- min_haz
      ud[i] <- Rx[i] * hd[i]
      next
    }
    
    S1_tgt <- lx[i + 1]
    U1_tgt <- lu_tgt[i + 1]
    
    solve_hd_given_hu <- function(hu_i) {
      f_hd <- function(hd_i) {
        P <- make_P(hu = hu_i, hd = hd_i, Rx_i = Rx[i], dt = age_int)
        nxt <- step_forward(H0, U0, P)
        (nxt["H"] + nxt["U"]) - S1_tgt
      }
      
      lo <- hd_bracket[1]
      hi <- hd_bracket[2]
      f_lo <- f_hd(lo)
      f_hi <- f_hd(hi)
      
      k <- 0
      while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 30) {
        hi <- hi * 2
        f_hi <- f_hd(hi)
        k <- k + 1
      }
      if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) return(NA_real_)
      
      uniroot(f_hd, lower = lo, upper = hi, tol = tol, maxiter = maxiter)$root
    }
    
    f_hu <- function(hu_i) {
      hd_i <- solve_hd_given_hu(hu_i)
      if (!is.finite(hd_i)) return(NA_real_)
      
      P <- make_P(hu = hu_i, hd = hd_i, Rx_i = Rx[i], dt = age_int)
      nxt <- step_forward(H0, U0, P)
      nxt["U"] - U1_tgt
    }
    
    lo <- hu_bracket[1]
    hi <- hu_bracket[2]
    f_lo <- f_hu(lo)
    f_hi <- f_hu(hi)
    
    k <- 0
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 30) {
      hi <- hi * 2
      f_hi <- f_hu(hi)
      k <- k + 1
    }
    
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      if (verbose) {
        warning(
          sprintf("Could not bracket hu root at age=%s (interval %s->%s). Using min_haz.",
                  age[i], age[i], age[i + 1]),
          call. = FALSE
        )
      }
      hu[i] <- min_haz
      hd[i] <- solve_hd_given_hu(hu[i])
      if (!is.finite(hd[i])) hd[i] <- min_haz
      ud[i] <- Rx[i] * hd[i]
      next
    }
    
    hu_i <- uniroot(f_hu, lower = lo, upper = hi, tol = tol, maxiter = maxiter)$root
    hd_i <- solve_hd_given_hu(hu_i)
    if (!is.finite(hd_i)) hd_i <- min_haz
    
    hu[i] <- pmax(hu_i, min_haz)
    hd[i] <- pmax(hd_i, min_haz)
    ud[i] <- pmax(Rx[i] * hd[i], min_haz)
  }
  
  extrap_last <- function(x) {
    if (is.finite(x[n - 1]) && is.finite(x[n - 2]) && x[n - 1] > 0 && x[n - 2] > 0) {
      exp(2 * log(x[n - 1]) - log(x[n - 2]))
    } else if (is.finite(x[n - 1])) {
      x[n - 1]
    } else {
      min_haz
    }
  }
  
  hu[n] <- pmax(extrap_last(hu), min_haz)
  hd[n] <- pmax(extrap_last(hd), min_haz)
  ud[n] <- pmax(extrap_last(ud), min_haz)
  uh[n] <- min_haz
  
  tibble::tibble(
    age = rep(age, times = 4),
    trans = rep(c("hu", "hd", "ud", "uh"), each = n),
    hazard = c(hu, hd, ud, uh)
  )
}