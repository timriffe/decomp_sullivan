
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
