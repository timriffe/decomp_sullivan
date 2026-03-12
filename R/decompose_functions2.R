# one-age converter -------------------------------------------------------

one_age_probs <- function(hu, hd, uh, ud, age_int = 1, system = c("returns", "noreturns")) {
  system <- match.arg(system)
  
  if (system == "noreturns") {
    lamH <- hu + hd
    one_minus_survH <- 1 - exp(-lamH * age_int)
    
    p_hu <- if (lamH > 0) (hu / lamH) * one_minus_survH else 0
    p_hd <- if (lamH > 0) (hd / lamH) * one_minus_survH else 0
    p_uh <- 0
    p_ud <- 1 - exp(-ud * age_int)
    
    return(c(hu = p_hu, hd = p_hd, uh = p_uh, ud = p_ud))
  }
  
  Q <- matrix(
    c(-(hu + hd),  hu,          hd,
      uh,        -(uh + ud),   ud,
      0,          0,           0),
    nrow = 3, byrow = TRUE
  )
  
  P <- expm::expm(Q * age_int)
  
  c(
    hu = P[1, 2],
    hd = P[1, 3],
    uh = P[2, 1],
    ud = P[2, 3]
  )
}


# full schedule -> 4 x A probability matrix -------------------------------

hazards_to_probmat <- function(hazard, trans, age, age_int = 1, system = c("returns", "noreturns")) {
  system <- match.arg(system)
  
  ages <- age[!duplicated(age)]
  A <- length(ages)
  
  haz_wide <- long_to_wide_transitions(age = age, trans = trans, hazard = hazard)
  
  out <- matrix(0, nrow = 4, ncol = A,
                dimnames = list(c("hu", "hd", "uh", "ud"), ages))
  
  for (a in seq_len(A)) {
    out[, a] <- one_age_probs(
      hu = haz_wide$hu[a],
      hd = haz_wide$hd[a],
      uh = haz_wide$uh[a],
      ud = haz_wide$ud[a],
      age_int = age_int,
      system = system
    )
  }
  
  out
}


# prob matrix -> transitions tibble ---------------------------------------

probmat_to_transitions <- function(probmat, ages) {
  tibble::tibble(
    age = rep(ages, each = 4L),
    from_to = rep(c("hu", "hd", "uh", "ud"), times = length(ages)),
    p = c(probmat)
  )
}


# prob matrix -> HLE ------------------------------------------------------

probmat_to_hle <- function(probmat, ages, age_int = 1, init) {
  transitions <- probmat_to_transitions(probmat = probmat, ages = ages)
  Lx <- mscalc:::calc_occupancies(
    transitions = transitions,
    age_interval = age_int,
    init = init,
    delim = ""
  )
  sum(Lx$h)
}


# tailored horiuchi for hazard -> HLE -------------------------------------

horiuchi_haz_cached <- function(pars1, pars2, trans, age, age_int = 1,
                                init, system = c("returns", "noreturns"),
                                N = 20) {
  system <- match.arg(system)
  
  stopifnot(length(pars1) == length(pars2))
  stopifnot(length(pars1) == length(trans))
  stopifnot(length(pars1) == length(age))
  
  d <- pars2 - pars1
  n <- length(pars1)
  delta <- d / N
  
  grad <- matrix(
    rep((0.5:(N - 0.5)) / N, n),
    byrow = TRUE,
    ncol = N
  )
  
  x <- pars1 + d * grad
  rownames(x) <- names(pars1)
  
  ages <- age[!duplicated(age)]
  A <- length(ages)
  
  # map each parameter index to:
  #   age position in 1:A
  #   transition slot in 1:4
  age_index <- match(age, ages)
  trans_index <- match(trans, c("hu", "hd", "uh", "ud"))
  
  cc <- matrix(0, nrow = n, ncol = N)
  
  for (j in seq_len(N)) {
    theta_j <- x[, j]
    
    # baseline probability schedule for this interpolation point
    prob0 <- hazards_to_probmat(
      hazard = theta_j,
      trans = trans,
      age = age,
      age_int = age_int,
      system = system
    )
    
    # wide hazards at this interpolation point, so we can edit one age block
    hw <- long_to_wide_transitions(age = age, trans = trans, hazard = theta_j)
    
    for (i in seq_len(n)) {
      a <- age_index[i]
      t <- trans_index[i]
      
      # grab current hazards for that age
      hu0 <- hw$hu[a]
      hd0 <- hw$hd[a]
      uh0 <- hw$uh[a]
      ud0 <- hw$ud[a]
      
      # +/- perturbation only for the touched hazard
      bump <- delta[i] / 2
      
      hplus  <- c(hu0, hd0, uh0, ud0)
      hminus <- c(hu0, hd0, uh0, ud0)
      
      hplus[t]  <- hplus[t]  + bump
      hminus[t] <- hminus[t] - bump
      
      # recompute only that one age block
      pplus_block <- one_age_probs(
        hu = hplus[1], hd = hplus[2], uh = hplus[3], ud = hplus[4],
        age_int = age_int, system = system
      )
      
      pminus_block <- one_age_probs(
        hu = hminus[1], hd = hminus[2], uh = hminus[3], ud = hminus[4],
        age_int = age_int, system = system
      )
      
      prob_plus <- prob0
      prob_minus <- prob0
      
      prob_plus[, a] <- pplus_block
      prob_minus[, a] <- pminus_block
      
      cc[i, j] <-
        probmat_to_hle(prob_plus, ages = ages, age_int = age_int, init = init) -
        probmat_to_hle(prob_minus, ages = ages, age_int = age_int, init = init)
    }
  }
  
  out <- rowSums(cc)
  names(out) <- names(pars1)
  out
}