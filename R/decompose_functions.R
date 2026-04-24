
long_to_wide_transitions <- function(age, trans, hazard) {
  ages <- age[!duplicated(age)]
  n <- length(ages)
  
  hu <- hd <- uh <- ud <- numeric(n)
  
  ai <- match(age, ages)
  
  # match trans once (no parsing)
  ti <- match(trans, c("hu","hd","uh","ud"))
  # scatter
  for (k in seq_along(hazard)) {
    i <- ai[k]
    t <- ti[k]
    if (t == 1L) hu[i] <- hazard[k]
    else if (t == 2L) hd[i] <- hazard[k]
    else if (t == 3L) uh[i] <- hazard[k]
    else if (t == 4L) ud[i] <- hazard[k]
    else stop("Unexpected trans value: ", trans[k])
  }
  
  list(age = ages, hu = hu, hd = hd, uh = uh, ud = ud)
}
haz_to_transitions_returns <- function(age, hu, hd, uh, ud, age_int = 1) {
  n <- length(age)
  stopifnot(length(hu) == n, length(hd) == n, length(uh) == n, length(ud) == n)
  
  out_age <- rep(age, each = 4L)
  out_tr  <- rep(c("hu","hd","uh","ud"), times = n)
  out_p   <- numeric(4L * n)
  
  for (i in seq_len(n)) {
    Q <- matrix(
      c(
        -(hu[i] + hd[i]),  hu[i],           hd[i],
        uh[i],           -(uh[i] + ud[i]), ud[i],
        0,                0,               0
      ),
      nrow = 3L, byrow = TRUE
    )
    
    P <- expm::expm(Q * age_int)
    
    k <- (i - 1L) * 4L
    out_p[k + 1L] <- P[1L, 2L] # hu
    out_p[k + 2L] <- P[1L, 3L] # hd
    out_p[k + 3L] <- P[2L, 1L] # uh
    out_p[k + 4L] <- P[2L, 3L] # ud
  }
  
  tibble::tibble(age = out_age, from_to = out_tr, p = out_p)
}

# identical result to above if uh=0, faster cuz no expm
haz_to_transitions_noreturns_cr <- function(age, hu, hd, ud, age_int = 1) {
  n <- length(age)
  stopifnot(length(hu) == n, length(hd) == n, length(ud) == n)
  
  lamH <- hu + hd
  # avoid 0/0 when lamH = 0
  one_minus_survH <- 1 - exp(-lamH * age_int)
  
  p_hu <- ifelse(lamH > 0, (hu / lamH) * one_minus_survH, 0)
  p_hd <- ifelse(lamH > 0, (hd / lamH) * one_minus_survH, 0)
  
  p_ud <- 1 - exp(-ud * age_int)
  p_uh <- rep(0, n)  # explicitly include as zeros for shape-compatibility
  
  tibble::tibble(
    age = rep(age, each = 4L),
    from_to = rep(c("hu","hd","uh","ud"), times = n),
    p = c(rbind(p_hu, p_hd, p_uh, p_ud))
  )
}

haz_to_hle <- function(hazard, trans, age, age_int = 1, init, system) {
  haz_wide <- long_to_wide_transitions(age, trans, hazard)
  if (system == "noreturns"){
    transitions <- haz_to_transitions_noreturns_cr(
      age = haz_wide$age, 
      hu = haz_wide$hu, 
      hd = haz_wide$hd, 
      ud = haz_wide$ud, 
      age_int = age_int
    )
  }
  if (system == "returns"){
    transitions <- haz_to_transitions_returns(
      age = haz_wide$age, 
      hu = haz_wide$hu, 
      hd = haz_wide$hd, 
      uh = haz_wide$uh, 
      ud = haz_wide$ud, 
      age_int = age_int
    )
  }
  Lx = mscalc:::calc_occupancies(transitions = transitions, age_interval = age_int, init = init, delim="")
  sum(Lx$h)
}



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


probmat_to_transitions <- function(probmat, ages) {
  tibble::tibble(
    age = rep(ages, each = 4L),
    from_to = rep(c("hu", "hd", "uh", "ud"), times = length(ages)),
    p = c(probmat)
  )
}

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
      ti <- trans_index[i]
      
      # grab current hazards for that age
      hu0 <- hw$hu[a]
      hd0 <- hw$hd[a]
      uh0 <- hw$uh[a]
      ud0 <- hw$ud[a]
      
      # +/- perturbation only for the touched hazard
      bump <- delta[i] / 2
      
      hplus  <- c(hu0, hd0, uh0, ud0)
      hminus <- c(hu0, hd0, uh0, ud0)
      
      hplus[ti]  <- hplus[ti]  + bump
      hminus[ti] <- hminus[ti] - bump
      
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