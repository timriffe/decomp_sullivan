
# --------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(expm)
  library(furrr)
  library(future)
  library(parallel)
})
validate_par_vec <- function(par_vec, piecewise = NULL) {
  if (!is.numeric(par_vec)) stop("par_vec must be numeric.", call. = FALSE)
  if (!(length(par_vec) %in% c(8, 12))) stop("par_vec must have length 8 (linear) or 12 (piecewise).", call. = FALSE)
  
  if (is.null(piecewise)) piecewise <- (length(par_vec) == 12)
  want <- ms_par_names(piecewise)
  
  # If unnamed, assign by position in the chosen representation
  if (is.null(names(par_vec)) || !all(want %in% names(par_vec))) {
    # allow legacy 8-vector even when piecewise=TRUE: we'll convert later
    if (length(par_vec) == 8) {
      names(par_vec) <- ms_par_names(FALSE)
    } else {
      names(par_vec) <- ms_par_names(TRUE)
    }
  }
  
  if (isTRUE(piecewise) && length(par_vec) == 12) {
    return(par_vec[ms_par_names(TRUE)])
  }
  if (!isTRUE(piecewise) && length(par_vec) == 8) {
    return(par_vec[ms_par_names(FALSE)])
  }
  
  # If mismatch (e.g., asked piecewise but got 8), 
  # just return canonical for its length.
  if (length(par_vec) == 8) return(par_vec[ms_par_names(FALSE)])
  par_vec[ms_par_names(TRUE)]
}

ms_par_names <- function(piecewise = FALSE) {
  if (!isTRUE(piecewise)) {
    c("i_hu","i_hd","i_uh","i_ud","s_hu","s_hd","s_uh","s_ud")
  } else {
    c(
      "i_hu","i_hd","i_uh","i_ud",
      "s1_hu","s1_hd","s1_uh","s1_ud", # first part of slope
      "s2_hu","s2_hd","s2_uh","s2_ud"  # second part of slope
    )
  }
}

as_par_vec <- function(pars) {
  if (is.data.frame(pars) || inherits(pars, "tbl_df")) {
    pars <- dplyr::as_tibble(pars)
    if (all(c("trans","intercept","slope") %in% names(pars))) {
      want <- c("hu","hd","uh","ud")
      pars <- dplyr::distinct(pars, trans, .keep_all = TRUE) |>
        dplyr::filter(.data$trans %in% want) |>
        dplyr::arrange(match(.data$trans, want))
      out <- c(
        setNames(as.numeric(pars$intercept), paste0("i_", pars$trans)),
        setNames(as.numeric(pars$slope),     paste0("s_", pars$trans))
      )
      return(validate_par_vec(out, piecewise = FALSE))
    }
    
    if (all(c("trans","intercept","slope1","slope2") %in% names(pars))) {
      want <- c("hu","hd","uh","ud")
      pars <- dplyr::distinct(pars, trans, .keep_all = TRUE) |>
        dplyr::filter(.data$trans %in% want) |>
        dplyr::arrange(match(.data$trans, want))
      out <- c(
        setNames(as.numeric(pars$intercept), paste0("i_", pars$trans)),
        setNames(as.numeric(pars$slope1),    paste0("s1_", pars$trans)),
        setNames(as.numeric(pars$slope2),    paste0("s2_", pars$trans))
      )
      return(validate_par_vec(out, piecewise = TRUE))
    }
    
    stop("pars df must have (trans, intercept, slope) or (trans, intercept, slope1, slope2).", call. = FALSE)
  }
  
  validate_par_vec(pars, piecewise = NULL)
}

as_piecewise_pars <- function(pars) {
  v <- validate_par_vec(as_par_vec(pars), piecewise = NULL)
  if (length(v) == 12) return(v)
  
  out <- c(
    v[c("i_hu","i_hd","i_uh","i_ud")],
    setNames(as.numeric(v[c("s_hu","s_hd","s_uh","s_ud")]), paste0("s1_", c("hu","hd","uh","ud"))),
    setNames(as.numeric(v[c("s_hu","s_hd","s_uh","s_ud")]), paste0("s2_", c("hu","hd","uh","ud")))
  )
  out[ms_par_names(TRUE)]
}

expand_hazards <- function(pars,
                           age = 50:100,
                           age_int = 1,
                           haz_shape = c("loglinear","piecewise","softkink"),
                           pivot_age = 75,
                           shape_width = 5,
                           knot_age = NULL) {
  
  # Backward compatibility: allow older code to pass knot_age
  if (!is.null(knot_age)) pivot_age <- knot_age
  
  haz_shape <- match.arg(haz_shape)
  
  v <- as_par_vec(pars)
  age_mid <- age + age_int / 2
  
  if (length(v) == 8) {
    # --- 8-param log-linear hazards (intercept + slope) ---
    pars_df <- pars_vec_to_df(v)
    
    tidyr::expand_grid(age = age, trans = pars_df$trans) |>
      dplyr::left_join(pars_df, by = "trans") |>
      dplyr::mutate(
        log_hazard = intercept + slope * (age + age_int/2 - min(age)),
        hazard = exp(log_hazard)
      ) |>
      dplyr::select(age, trans, hazard)
    
  } else {
    # --- 12-param piecewise family: i_*, s1_*, s2_* ---
    i  <- v[c("i_hu","i_hd","i_uh","i_ud")]
    s1 <- v[c("s1_hu","s1_hd","s1_uh","s1_ud")]
    s2 <- v[c("s2_hu","s2_hd","s2_uh","s2_ud")]
    
    trans <- c("hu","hd","uh","ud")
    i  <- setNames(as.numeric(i),  trans)
    s1 <- setNames(as.numeric(s1), trans)
    s2 <- setNames(as.numeric(s2), trans)
    
    # centered age for stability
    a0 <- min(age)
    x  <- (age_mid - a0)
    
    # hard-kink term (piecewise)
    kink_hard <- pmax(0, age_mid - pivot_age)
    
    # smooth-kink term: softplus((a - pivot)/w), anchored at first interval mid-age
    # This yields a C^1 curve (no sharp kink), and (s2-s1)=0 collapses to log-linear.
    softplus <- function(z) log1p(exp(z))
    z  <- (age_mid - pivot_age) / shape_width
    z0 <- (min(age_mid) - pivot_age) / shape_width
    kink_soft <- shape_width * (softplus(z) - softplus(z0))
    
    tidyr::expand_grid(age = age, trans = trans) |>
      dplyr::mutate(
        age_mid   = age + age_int/2,
        x         = age_mid - a0,
        kink_term =
          if (haz_shape == "piecewise") {
            pmax(0, age_mid - pivot_age)
          } else if (haz_shape == "softkink") {
            shape_width * (softplus((age_mid - pivot_age)/shape_width) - softplus(z0))
          } else {
            0
          },
        intercept = i[trans],
        slope1    = s1[trans],
        slope2    = s2[trans],
        log_hazard = intercept + slope1 * x + (slope2 - slope1) * kink_term,
        hazard     = exp(log_hazard)
      ) |>
      dplyr::select(age, trans, hazard)
  }
}

run_mslt <- function(pars,
                     age,
                     init = c(H = 1, U = 0),
                     age_int = 1,
                     haz_shape = c("loglinear","piecewise","softkink"),
                     pivot_age = 75,
                     shape_width = 5,
                     knot_age = NULL) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  haz <- expand_hazards(pars, age = age, age_int = age_int, haz_shape = haz_shape,
                        pivot_age = pivot_age, shape_width = shape_width)
  P   <- haz_to_probs(haz, age = age, age_int = age_int)
  calculate_lt(P, init = init, age_int = age_int)
}

mslt_summary_vector <- function(pars,
                                age,
                                init = c(H = 1, U = 0),
                                age_int = 1,
                                haz_shape = c("loglinear","piecewise","softkink"),
                                pivot_age = 75,
                                shape_width = 5,
                                knot_age = NULL,
                                outputs = c("lx", "prevalence"),
                                weights = NULL,
                                standardize = TRUE) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  
  
  lt <- run_mslt(pars, age = age, init = init, age_int = age_int, haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width)
  
  if (!all(outputs %in% names(lt))) {
    stop("Some requested outputs are not in lifetable: ",
         paste(setdiff(outputs, names(lt)), collapse = ", "),
         call. = FALSE)
  }
  
  blocks <- lapply(outputs, function(nm) as.numeric(lt[[nm]]))
  
  if (!is.null(weights)) {
    if (!is.list(weights) || is.null(names(weights))) {
      stop("weights must be a named list keyed by output names.", call. = FALSE)
    }
    for (nm in outputs) {
      if (!is.null(weights[[nm]])) {
        w <- weights[[nm]]
        if (length(w) == 1) w <- rep(w, length(age))
        if (length(w) != length(age)) {
          stop("weights[['", nm, "']] must be length 1 or length(age).", call. = FALSE)
        }
        blocks[[which(outputs == nm)]] <- blocks[[which(outputs == nm)]] * as.numeric(w)
      }
    }
  }
  
  if (standardize) {
    blocks <- lapply(blocks, function(x) {
      s <- stats::sd(x)
      if (is.finite(s) && s > 0) (x - mean(x)) / s else x - mean(x)
    })
  }
  
  unlist(blocks, use.names = FALSE)
}

mslt_jacobian <- function(base_pars,
                          free_names,
                          age,
                          init = c(H = 1, U = 0),
                          age_int = 1,
                          haz_shape = c("loglinear","piecewise","softkink"),
                          pivot_age = 75,
                          shape_width = 5,
                          knot_age = NULL,
                          outputs = c("lx", "prevalence"),
                          weights = NULL,
                          standardize = TRUE,
                          eps = 1e-4,
                          method = c("central", "forward"),
                          step_scale = c("relative", "absolute")) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  
  
  method     <- match.arg(method)
  step_scale <- match.arg(step_scale)
  
  theta0 <- as_par_vec(base_pars)
  
  y0 <- mslt_summary_vector(theta0, age, init, age_int, haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width, knot_age = knot_age, outputs, weights, standardize)
  m  <- length(y0)
  p  <- length(free_names)
  
  J <- matrix(NA_real_, nrow = m, ncol = p, dimnames = list(NULL, free_names))
  
  for (j in seq_len(p)) {
    nm <- free_names[j]
    
    step <- if (step_scale == "relative") eps * (abs(theta0[[nm]]) + 1) else eps
    
    theta_plus <- theta0
    theta_plus[[nm]] <- theta0[[nm]] + step
    y_plus <- mslt_summary_vector(theta_plus, age, init, age_int, haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width, knot_age = knot_age, outputs, weights, standardize)
    
    if (method == "central") {
      theta_minus <- theta0
      theta_minus[[nm]] <- theta0[[nm]] - step
      y_minus <- mslt_summary_vector(theta_minus, age, init, age_int, haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width, knot_age = knot_age, outputs, weights, standardize)
      J[, j] <- (y_plus - y_minus) / (2 * step)
    } else {
      J[, j] <- (y_plus - y0) / step
    }
  }
  
  list(J = J, y0 = y0, theta0 = theta0[free_names], free_names = free_names,
       meta = list(age = age, age_int = age_int, haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width, outputs = outputs,
                   standardize=standardize, eps=eps, method=method, step_scale=step_scale))
}

null_directions_from_jacobian <- function(J, k = NULL, tol = 1e-3) {
  sv <- svd(J)
  s  <- sv$d
  V  <- sv$v
  
  if (!is.null(k)) {
    k <- max(1, min(k, ncol(J)))
    idx <- (length(s) - k + 1):length(s)
  } else {
    cutoff <- tol * max(s)
    idx <- which(s <= cutoff)
    if (length(idx) == 0) idx <- length(s)
  }
  
  idx <- idx[order(s[idx])] # most null first
  
  dirs <- V[, idx, drop = FALSE]
  rownames(dirs) <- colnames(J)
  colnames(dirs) <- paste0("dir", seq_len(ncol(dirs)))
  
  list(directions = dirs, singular_values = s, svd = sv)
}

find_null_directions <- function(base_pars,
                                 free_names,
                                 age,
                                 init = c(H = 1, U = 0),
                                 age_int = 1,
                                 haz_shape = c("loglinear","piecewise","softkink"),
                                 pivot_age = 75,
                                 shape_width = 5,
                                 knot_age = NULL,
                                 outputs = c("lx","prevalence"),
                                 tol = 1e-3,
                                 eps = 1e-4,
                                 method = c("central","forward"),
                                 step_scale = c("relative","absolute"),
                                 standardize = TRUE) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  
  
  jac <- mslt_jacobian(
    base_pars = base_pars,
    free_names = free_names,
    age = age,
    init = init,
    age_int = age_int,
    haz_shape = haz_shape,
    pivot_age = pivot_age,
    shape_width = shape_width,
    knot_age = knot_age,
    outputs = outputs,
    standardize = standardize,
    eps = eps,
    method = match.arg(method),
    step_scale = match.arg(step_scale)
  )
  
  nd <- null_directions_from_jacobian(jac$J, tol = tol)
  
  list(
    base_pars = as_par_vec(base_pars),
    free_names = free_names,
    jacobian = jac,
    directions = nd$directions,
    singular_values = nd$singular_values,
    svd = nd$svd,
    meta = jac$meta
  )
}

apply_direction <- function(base_pars, v, step, free_names = NULL) {
  theta <- as_par_vec(base_pars)
  
  v <- as.numeric(v)
  
  if (is.null(free_names)) free_names <- names(v)
  if (is.null(free_names) || any(free_names == "")) {
    stop("apply_direction(): provide free_names (names(v) missing).", call. = FALSE)
  }
  
  if (length(v) != length(free_names)) stop("Length mismatch: length(v) != length(free_names).", call. = FALSE)
  
  theta[free_names] <- theta[free_names] + step * v
  theta
}

apply_directions <- function(base_pars, V, steps, free_names = NULL) {
  theta <- as_par_vec(base_pars)
  V <- as.matrix(V)
  steps <- as.numeric(steps)
  
  if (length(steps) != ncol(V)) stop("steps must have length ncol(V).", call. = FALSE)
  
  if (is.null(free_names)) free_names <- rownames(V)
  if (is.null(free_names)) stop("Need free_names or rownames(V).", call. = FALSE)
  
  if (nrow(V) != length(free_names)) stop("nrow(V) must match length(free_names).", call. = FALSE)
  
  theta[free_names] <- theta[free_names] + as.numeric(V %*% steps)
  theta
}

make_bounds <- function(free_names,
                        slope_floor = 1e-6,
                        slope_ceiling = Inf) {
  
  lo <- rep(-Inf, length(free_names)); hi <- rep( Inf, length(free_names))
  names(lo) <- names(hi) <- free_names
  
  # linear slopes
  if ("s_hu" %in% free_names) { lo["s_hu"] <- slope_floor; hi["s_hu"] <- slope_ceiling }
  if ("s_hd" %in% free_names) { lo["s_hd"] <- slope_floor; hi["s_hd"] <- slope_ceiling }
  if ("s_ud" %in% free_names) { lo["s_ud"] <- slope_floor; hi["s_ud"] <- slope_ceiling }
  if ("s_uh" %in% free_names) { lo["s_uh"] <- -slope_ceiling; hi["s_uh"] <- -slope_floor }
  
  # piecewise slopes
  for (nm in c("s1_hu","s1_hd","s1_ud","s2_hu","s2_hd","s2_ud")) {
    if (nm %in% free_names) { lo[nm] <- slope_floor; hi[nm] <- slope_ceiling }
  }
  for (nm in c("s1_uh","s2_uh")) {
    if (nm %in% free_names) { lo[nm] <- -slope_ceiling; hi[nm] <- -slope_floor }
  }
  
  list(lower = unname(lo), upper = unname(hi))
}

mslt_objective_ridge <- function(par_free,
                                 free_names,
                                 base_pars,
                                 ground_summary,
                                 age,
                                 init = c(H=1,U=0),
                                 age_int = 1,
                                 haz_shape = c("loglinear","piecewise","softkink"),
                                 pivot_age = 75,
                                 shape_width = 5,
                                 knot_age = NULL,
                                 fixed = NULL,
                                 w_lx = 1,
                                 w_prev = 1,
                                 lambda = 0.3,
                                 theta_target = NULL,
                                 par_scale = NULL,
                                 lambda_kink = 0,
                                 kink_trans = c("hu","hd","uh","ud"),
                                 # ---- optional: separate age-weights ----
                                 age_weight_lx = c("none","lx","Lx"),
                                 age_weight_prev = c("none","lx","Lx"),
                                 prev_band = c(52, 70),
                                 prev_band_mult = 1,
                                 # ---- NEW: optional integrand loss ----
                                 w_int = 0,
                                 integrand = c("health","unhealthy"),
                                 age_weight_int = c("none","lx","Lx")) {
  
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  age_weight_lx   <- match.arg(age_weight_lx)
  age_weight_prev <- match.arg(age_weight_prev)
  integrand <- match.arg(integrand)
  age_weight_int <- match.arg(age_weight_int)
  
  if (length(par_free) != length(free_names)) stop("par_free and free_names must have same length.", call. = FALSE)
  
  cand <- as_par_vec(base_pars)
  
  if (!is.null(fixed)) {
    if (is.null(names(fixed)) || any(names(fixed) == "")) stop("fixed must be named.", call. = FALSE)
    cand[names(fixed)] <- as.numeric(fixed)
  }
  cand[free_names] <- as.numeric(par_free)
  
  lt_mod <- run_mslt(
    cand, age = age, init = init, age_int = age_int,
    haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width
  ) |>
    dplyr::select(age, lx1 = lx, prevalence1 = prevalence) |>
    dplyr::left_join(ground_summary, by = "age") |>
    dplyr::mutate(
      loss_lx   = (lx1 - lx)^2,
      loss_prev = (prevalence1 - prevalence)^2
    )
  
  # ---- helper: build age weights from ground curves ----
  build_w <- function(which) {
    if (which == "none") return(rep(1, nrow(lt_mod)))
    if (which == "lx")   return(pmax(1e-12, lt_mod$lx))
    # "Lx" proxy from *ground* lx
    lxg <- pmax(1e-12, lt_mod$lx)
    lxg_next <- dplyr::lead(lxg, default = dplyr::last(lxg))
    pmax(1e-12, (lxg + lxg_next) / 2 * age_int)
  }
  
  w_age_lx   <- build_w(age_weight_lx)
  w_age_prev <- build_w(age_weight_prev)
  w_age_int  <- build_w(age_weight_int)
  
  # ---- prevalence band multiplier ----
  w_prev_age <- rep(1, nrow(lt_mod))
  if (!is.null(prev_band) && length(prev_band) == 2 && is.finite(prev_band_mult) && prev_band_mult > 1) {
    lo <- min(prev_band); hi <- max(prev_band)
    w_prev_age[lt_mod$age >= lo & lt_mod$age <= hi] <- prev_band_mult
  }
  
  # ---- optional integrand term: (1-prev)*Lx or prev*Lx ----
  int_loss <- 0
  if (!is.null(w_int) && is.finite(w_int) && w_int > 0) {
    # Lx proxies from model and ground
    lx1  <- pmax(1e-12, lt_mod$lx1)
    lx1n <- dplyr::lead(lx1, default = dplyr::last(lx1))
    Lx1  <- (lx1 + lx1n) / 2 * age_int
    
    lxg  <- pmax(1e-12, lt_mod$lx)
    lxgn <- dplyr::lead(lxg, default = dplyr::last(lxg))
    Lxg  <- (lxg + lxgn) / 2 * age_int
    
    if (integrand == "health") {
      I1 <- (1 - lt_mod$prevalence1) * Lx1
      Ig <- (1 - lt_mod$prevalence)  * Lxg
    } else {
      I1 <- (lt_mod$prevalence1) * Lx1
      Ig <- (lt_mod$prevalence)  * Lxg
    }
    
    lt_mod$loss_int <- (I1 - Ig)^2
    int_loss <- sum(w_age_int * lt_mod$loss_int, na.rm = TRUE)
  }
  
  fit_loss <- sum(
    w_lx   * w_age_lx   * lt_mod$loss_lx +
      w_prev * w_age_prev * w_prev_age * lt_mod$loss_prev,
    na.rm = TRUE
  ) + w_int * int_loss
  
  # ---- ridge identity term ----
  reg_loss <- 0
  if (!is.null(lambda) && lambda > 0) {
    if (is.null(theta_target)) theta_target <- as_par_vec(base_pars)
    theta_target <- as_par_vec(theta_target)
    
    if (is.null(par_scale)) {
      par_scale <- rep(1, length(free_names))
      names(par_scale) <- free_names
      for (nm in free_names) par_scale[nm] <- max(1e-6, abs(theta_target[[nm]]) + 1)
    } else {
      if (is.null(names(par_scale))) names(par_scale) <- free_names
      par_scale <- par_scale[free_names]
      par_scale[!is.finite(par_scale) | par_scale <= 0] <- 1
    }
    
    diff <- (cand[free_names] - theta_target[free_names]) / par_scale
    reg_loss <- lambda * sum(as.numeric(diff)^2)
  }
  
  # ---- kink penalty: only for kink_trans ----
  kink_loss <- 0
  if (!is.null(lambda_kink) && lambda_kink > 0 && length(kink_trans) > 0) {
    cand_pw <- if (length(cand) == 12) cand else as_piecewise_pars(cand)
    for (tr in kink_trans) {
      s1n <- paste0("s1_", tr)
      s2n <- paste0("s2_", tr)
      if (s1n %in% names(cand_pw) && s2n %in% names(cand_pw)) {
        kink_loss <- kink_loss + (cand_pw[[s2n]] - cand_pw[[s1n]])^2
      }
    }
    kink_loss <- lambda_kink * kink_loss
  }
  
  fit_loss + reg_loss + kink_loss
}

polish_candidate_ridge <- function(theta_start,
                                   ground_summary,
                                   free_polish,
                                   age,
                                   init = c(H=1,U=0),
                                   age_int = 1,
                                   haz_shape = c("loglinear","piecewise","softkink"),
                                   pivot_age = 75,
                                   shape_width = 5,
                                   knot_age = NULL,
                                   fixed = NULL,
                                   lambda = 0.3,
                                   theta_target = theta_start,
                                   slope_floor = 1e-6,
                                   maxit1 = 1500,
                                   lambda2 = NULL,
                                   maxit2 = 800,
                                   # kink, pass-specific
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
                                   # optional pass-2 convergence tightening
                                   factr1 = 1e6,
                                   factr2 = 1e6) {
  
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  age_weight_lx <- match.arg(age_weight_lx)
  age_weight_prev <- match.arg(age_weight_prev)
  integrand <- match.arg(integrand)
  age_weight_int <- match.arg(age_weight_int)
  
  b <- make_bounds(free_polish, slope_floor = slope_floor)
  
  run_one <- function(theta_start_local, theta_target_local, lambda_local, maxit_local,
                      lambda_kink_local, kink_trans_local, factr_local) {
    
    res <- optim(
      par = as.numeric(theta_start_local[free_polish]),
      fn  = mslt_objective_ridge,
      method = "L-BFGS-B",
      lower = b$lower,
      upper = b$upper,
      control = list(maxit = maxit_local, factr = factr_local),
      free_names = free_polish,
      base_pars = theta_start_local,
      ground_summary = ground_summary,
      age = age,
      init = init,
      age_int = age_int,
      haz_shape = haz_shape,
      pivot_age = pivot_age,
      shape_width = shape_width,
      knot_age = knot_age,
      fixed = fixed,
      lambda = lambda_local,
      theta_target = theta_target_local,
      lambda_kink = lambda_kink_local,
      kink_trans = kink_trans_local,
      w_lx = w_lx,
      w_prev = w_prev,
      age_weight_lx = age_weight_lx,
      age_weight_prev = age_weight_prev,
      prev_band = prev_band,
      prev_band_mult = prev_band_mult,
      w_int = w_int,
      integrand = integrand,
      age_weight_int = age_weight_int
    )
    
    cand <- as_par_vec(theta_start_local)
    if (!is.null(fixed)) cand[names(fixed)] <- as.numeric(fixed)
    cand[free_polish] <- res$par
    cand
  }
  
  theta1 <- run_one(theta_start, theta_target, lambda, maxit1,
                    lambda_kink1, kink_trans1, factr1)
  
  if (!is.null(lambda2) && is.finite(lambda2) && lambda2 >= 0) {
    theta2 <- run_one(theta1, theta1, lambda2, maxit2,
                      lambda_kink2, kink_trans2, factr2)
    return(theta2)
  }
  
  theta1
}

polish_many <- function(theta0_list,
                        ground_summary,
                        free_polish,
                        age,
                        init = c(H=1,U=0),
                        age_int = 1,
                        haz_shape = c("loglinear","piecewise","softkink"),
                        pivot_age = 75,
                        shape_width = 5,
                        knot_age = NULL,
                        fixed = NULL,
                        # stage 1 ridge
                        lambda = 0.3,
                        lambda_kink1 = 0,
                        kink_trans1 = c("hu","hd","uh","ud"),
                        maxit1 = 800,
                        factr1 = 1e6,
                        # optional stage 2 ridge
                        lambda2 = NULL,
                        lambda_kink2 = 0,
                        kink_trans2 = c("hu","hd","uh","ud"),
                        maxit2 = 1500,
                        factr2 = 1e4,
                        # fit weights / weighting schemes
                        w_lx = 1,
                        w_prev = 1,
                        age_weight_lx = c("none","lx","Lx"),
                        age_weight_prev = c("none","lx","Lx"),
                        prev_band = NULL,
                        prev_band_mult = 1,
                        # parallel
                        n_cores = 1) {
  
  haz_shape <- match.arg(haz_shape)
  age_weight_lx   <- match.arg(age_weight_lx)
  age_weight_prev <- match.arg(age_weight_prev)
  
  one <- function(th0) {
    polish_candidate_ridge(
      theta_start = th0,
      theta_target = th0,           # anchor to candidate
      ground_summary = ground_summary,
      free_polish = free_polish,
      age = age,
      init = init,
      age_int = age_int,
      haz_shape = haz_shape,
      pivot_age = pivot_age,
      shape_width = shape_width,
      knot_age = knot_age,
      fixed = fixed,
      # stage 1
      lambda = lambda,
      lambda_kink1 = lambda_kink1,
      kink_trans1 = kink_trans1,
      maxit1 = maxit1,
      factr1 = factr1,
      # stage 2 (optional)
      lambda2 = lambda2,
      lambda_kink2 = lambda_kink2,
      kink_trans2 = kink_trans2,
      maxit2 = maxit2,
      factr2 = factr2,
      # fit weights
      w_lx = w_lx,
      w_prev = w_prev,
      age_weight_lx = age_weight_lx,
      age_weight_prev = age_weight_prev,
      prev_band = prev_band,
      prev_band_mult = prev_band_mult
    )
  }
  
  if (is.null(n_cores) || n_cores <= 1) {
    return(lapply(theta0_list, one))
  }
  
  # parallel
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Install packages 'future' and 'future.apply' for n_cores > 1.", call. = FALSE)
  }
  
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  
  future::plan(future::multisession, workers = n_cores)
  
  future.apply::future_lapply(theta0_list, one, future.seed = TRUE)
}

score_fit <- function(theta, ground_summary, age,
                      init = c(H=1,U=0),
                      age_int = 1,
                      haz_shape = c("loglinear","piecewise","softkink"),
                      pivot_age = 75,
                      shape_width = 5,
                      knot_age = NULL) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  lt <- run_mslt(theta, age = age, init = init, age_int = age_int,
                 haz_shape = haz_shape, pivot_age = pivot_age, shape_width = shape_width)
  lt |>
    dplyr::select(age, lx1=lx, prev1=prevalence) |>
    dplyr::left_join(ground_summary, by="age") |>
    dplyr::summarise(
      rmse_lx = sqrt(mean((lx1 - lx)^2)),
      rmse_prev = sqrt(mean((prev1 - prevalence)^2)),
      maxabs_lx = max(abs(lx1 - lx)),
      maxabs_prev = max(abs(prev1 - prevalence))
    )
}

make_Q <- function(x) {
  with(as.list(x),
       matrix(c(
         -(hu + hd),  hu,  hd,
         uh, -(uh + ud), ud,
         0,   0,   0
       ), nrow = 3, byrow = TRUE)
  )
}

haz_to_probs <- function(hazard, age, age_int = 1) {
  
  hazard <- hazard %>%
    dplyr::select(age, trans, hazard)
  
  out <- hazard %>%
    dplyr::group_nest(age) %>%
    dplyr::mutate(
      Q = purrr::map(.data$data, ~ make_Q(tibble::deframe(dplyr::select(.x, trans, hazard)))),
      P = purrr::map(.data$Q, ~ (.x * age_int) %>%
                       expm::expm() %>%
                       as.data.frame() %>%
                       rlang::set_names(c("H", "U", "D")) %>%
                       cbind("from" = c("H", "U", "D")))
    ) %>%
    dplyr::select(age, P) %>%
    tidyr::unnest(P)
  
  out
}
calculate_lt <- function(P, init = c(H = 1, U = 0), age_int = 1,
                         terminal = c( "carry","zero", "point")) {
  terminal <- match.arg(terminal)
  
  out <- P |>
    tidyr::pivot_longer(-c(age, from), names_to = "to", values_to = "p") |>
    dplyr::filter(.data$from != "D") |>
    tidyr::unite("state", c("from", "to"), sep = "") |>
    tidyr::pivot_wider(names_from = state, values_from = p) |>
    dplyr::group_modify(~{
      .x$lh <- numeric(nrow(.x))
      .x$lu <- numeric(nrow(.x))
      
      .x$lh[1] <- init["H"]
      .x$lu[1] <- init["U"]
      
      for (ii in seq_len(nrow(.x) - 1)) {
        .x$lh[ii + 1] <- .x$lh[ii] * .x$HH[ii] + .x$lu[ii] * .x$UH[ii]
        .x$lu[ii + 1] <- .x$lu[ii] * .x$UU[ii] + .x$lh[ii] * .x$HU[ii]
      }
      .x
    }) |>
    dplyr::mutate(
      lx = .data$lh + .data$lu,
      
      # point prevalence at exact age x
      prevalence_point = dplyr::if_else(.data$lx > 0, .data$lu / .data$lx, NA_real_),
      
      # interval person-years
      Lh = age_int * (.data$lh + dplyr::lead(.data$lh, default = 0)) / 2,
      Lu = age_int * (.data$lu + dplyr::lead(.data$lu, default = 0)) / 2,
      Lx = .data$Lh + .data$Lu
    )
  
  # interval prevalence, with explicit terminal handling
  prev_int <- with(out, Lu / Lx)
  
  n <- nrow(out)
  if (n >= 1) {
    if (terminal == "zero") {
      prev_int[n] <- NA_real_
    } else if (terminal == "point") {
      prev_int[n] <- out$prevalence_point[n]
    } else if (terminal == "carry") {
      prev_int[n] <- if (n >= 2) prev_int[n - 1] else out$prevalence_point[n]
    }
  }
  
  out |>
    dplyr::mutate(
      prevalence_interval = prev_int,
      # make interval prevalence the default prevalence
      prevalence = .data$prevalence_interval
    )
}


`%||%` <- function(a, b) if (!is.null(a)) a else b
clamp <- function(x, lo, hi) pmin(hi, pmax(lo, x))

interval_prev_trapezoid <- function(lx, prev, last = c("carry", "point")) {
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
  den <- lx[1:(n - 1)] + lx[2:n]
  out[1:(n - 1)] <- ifelse(den > 0, (lu[1:(n - 1)] + lu[2:n]) / den, NA_real_)
  out[n] <- if (last == "point") prev[n] else out[n - 1]
  out
}

.coerce_state_prev <- function(prev,
                               prevalence_point = NULL,
                               lx = NULL,
                               n = NULL) {
  out <- if (!is.null(prevalence_point)) prevalence_point else prev
  out <- as.numeric(out)
  
  if (!is.null(n) && length(out) != n) {
    stop("state prevalence vector must have length n.", call. = FALSE)
  }
  if (!is.null(lx) && length(out) != length(lx)) {
    stop("state prevalence vector must match length(lx).", call. = FALSE)
  }
  out
}

build_states_targets <- function(age, lx, prevalence_point) {
  age <- as.numeric(age)
  lx  <- as.numeric(lx)
  prev_pt <- as.numeric(prevalence_point)
  n <- length(age)
  
  if (length(lx) != n) {
    stop("lx must have same length as age.", call. = FALSE)
  }
  if (length(prev_pt) != n) {
    stop("prevalence_point must have same length as age.", call. = FALSE)
  }
  
  prev_pt <- pmin(pmax(prev_pt, 0), 1)
  
  # exact-age state stocks
  lH <- lx * (1 - prev_pt)
  lU <- lx * prev_pt
  lD <- 1 - lx
  
  # targets for age x+1 implied by the supplied exact-age stocks
  nextHU <- cbind(lH, lU)
  if (n >= 2) {
    nextHU[1:(n - 1), ] <- cbind(lH[2:n], lU[2:n])
    nextHU[n, ] <- nextHU[n - 1, ]
  }
  
  list(
    age = age,
    lx = lx,
    prev = prev_pt,
    lH = lH,
    lU = lU,
    lD = lD,
    nextHU = nextHU
  )
}
make_Q_returns <- function(hu, hd, uh, ud) {
  matrix(c(
    -(hu + hd),  hu,         hd,
    uh,        -(uh + ud),   ud,
    0,           0,          0
  ), nrow = 3, byrow = TRUE)
}

derive_noreturns_hazards <- function(
    age,
    lx,
    prev_pt,
    Rx,
    age_int = 1,
    min_haz = 1e-12,
    tol = 1e-14,
    maxiter = 100,
    hu_bracket = c(1e-12, 5),
    hd_bracket = c(1e-12, 5),
    verbose = TRUE
) {
  
  age <- as.numeric(age)
  lx  <- as.numeric(lx)
  prev_pt <- as.numeric(prev_pt)
  
  n <- length(age)
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  if (length(lx) != n) stop("lx must match age.", call. = FALSE)
  if (length(prev_pt) != n) stop("prev_pt must match age.", call. = FALSE)
  
  if (length(Rx) == 1) Rx <- rep(Rx, n)
  if (length(Rx) != n) stop("Rx must be length 1 or n.", call. = FALSE)
  
  Rx <- pmax(as.numeric(Rx), min_haz)
  
  # --- targets ---
  U <- lx * prev_pt
  H <- lx - U
  
  # --- helpers ---
  make_P <- function(hu, hd, ud, dt, min_haz) {
    uh <- min_haz
    Q <- matrix(c(
      -(hu + hd),  hu,  hd,
      uh, -(uh + ud), ud,
      0, 0, 0
    ), nrow = 3, byrow = TRUE)
    expm::expm(Q * dt)
  }
  
  step <- function(H0, U0, P) {
    v <- c(H0, U0, 0) %*% P
    c(H = as.numeric(v[1]), U = as.numeric(v[2]))
  }
  
  # --- storage ---
  hu <- numeric(n)
  hd <- numeric(n)
  ud <- numeric(n)
  uh <- rep(min_haz, n)
  
  # --- main loop ---
  for (i in seq_len(n - 1)) {
    
    H0 <- H[i]
    U0 <- U[i]
    
    if ((H0 + U0) <= min_haz) {
      hu[i] <- min_haz
      hd[i] <- min_haz
      ud[i] <- Rx[i] * hd[i]
      next
    }
    
    S1_tgt <- lx[i + 1]
    U1_tgt <- U[i + 1]
    
    # --- inner: solve hd for survival ---
    solve_hd <- function(hu_i) {
      
      f_hd <- function(hd_i) {
        ud_i <- Rx[i] * hd_i
        P <- make_P(hu_i, hd_i, ud_i, age_int, min_haz)
        nxt <- step(H0, U0, P)
        (nxt["H"] + nxt["U"]) - S1_tgt
      }
      
      lo <- hd_bracket[1]
      hi <- hd_bracket[2]
      
      f_lo <- f_hd(lo)
      f_hi <- f_hd(hi)
      
      k <- 0
      while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 25) {
        hi <- hi * 2
        f_hi <- f_hd(hi)
        k <- k + 1
      }
      
      if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
        return(NA_real_)
      }
      
      uniroot(f_hd, c(lo, hi), tol = tol, maxiter = maxiter)$root
    }
    
    # --- outer: solve hu for U ---
    f_hu <- function(hu_i) {
      hd_i <- solve_hd(hu_i)
      if (!is.finite(hd_i)) return(NA_real_)
      
      ud_i <- Rx[i] * hd_i
      P <- make_P(hu_i, hd_i, ud_i, age_int, min_haz)
      nxt <- step(H0, U0, P)
      nxt["U"] - U1_tgt
    }
    
    lo <- hu_bracket[1]
    hi <- hu_bracket[2]
    
    f_lo <- f_hu(lo)
    f_hi <- f_hu(hi)
    
    k <- 0
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 25) {
      hi <- hi * 2
      f_hi <- f_hu(hi)
      k <- k + 1
    }
    
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      if (verbose) {
        warning(sprintf("No root for hu at age %s", age[i]), call. = FALSE)
      }
      hu[i] <- min_haz
      hd[i] <- min_haz
      ud[i] <- Rx[i] * hd[i]
      next
    }
    
    hu_i <- uniroot(f_hu, c(lo, hi), tol = tol, maxiter = maxiter)$root
    hd_i <- solve_hd(hu_i)
    
    hu[i] <- pmax(hu_i, min_haz)
    hd[i] <- pmax(hd_i, min_haz)
    ud[i] <- pmax(Rx[i] * hd[i], min_haz)
  }
  
  # --- last age extrapolation ---
  extrap <- function(x) {
    if (is.finite(x[n-1]) && is.finite(x[n-2]) && x[n-1] > 0 && x[n-2] > 0) {
      exp(2 * log(x[n-1]) - log(x[n-2]))
    } else {
      x[n-1]
    }
  }
  
  hu[n] <- pmax(extrap(hu), min_haz)
  hd[n] <- pmax(extrap(hd), min_haz)
  ud[n] <- pmax(extrap(ud), min_haz)
  uh[n] <- min_haz
  
  tibble::tibble(
    age = rep(age, 4),
    trans = rep(c("hu", "hd", "ud", "uh"), each = n),
    hazard = c(hu, hd, ud, uh)
  )
}
make_candidates <- function(base_pars,
                            free_names,
                            V2,
                            grid_steps = c(-0.2, 0.2),
                            exclude_zero = TRUE) {
  base_pars <- as_par_vec(base_pars)
  V2 <- as.matrix(V2)
  stopifnot(ncol(V2) >= 2)
  if (is.null(rownames(V2))) rownames(V2) <- free_names
  if (is.null(free_names)) free_names <- rownames(V2)
  
  combos <- tidyr::expand_grid(s1 = grid_steps, s2 = grid_steps) |>
    dplyr::mutate(is_zero = (.data$s1 == 0 & .data$s2 == 0))
  if (isTRUE(exclude_zero)) combos <- combos |> dplyr::filter(!.data$is_zero)
  
  theta0 <- lapply(seq_len(nrow(combos)), function(i) {
    steps <- c(combos$s1[i], combos$s2[i])
    apply_directions(base_pars, V2[,1:2, drop=FALSE], steps = steps, free_names = free_names)
  })
  
  list(grid = combos, theta0 = theta0)
}



derive_returns_hazard_starts <- function(
    age,
    lx,
    prev_pt,
    Rx,
    theta_start,
    age_int = 1,
    pivot_age = 75,
    shape_width = 5,
    min_haz = 1e-12,
    maxit = 1500,
    factr = 1e6,
    w_lx = 1,
    w_prev = 1,
    slope_floor = 1e-6,
    verbose = TRUE
) {
  
  age <- as.numeric(age)
  lx <- as.numeric(lx)
  prev_pt <- as.numeric(prev_pt)
  n <- length(age)
  
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  if (length(lx) != n) stop("lx must match age.", call. = FALSE)
  if (length(prev_pt) != n) stop("prev_pt must match age.", call. = FALSE)
  
  if (length(Rx) == 1) Rx <- rep(as.numeric(Rx), n)
  Rx <- as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or n.", call. = FALSE)
  Rx <- pmax(Rx, min_haz)
  
  # ---- exact no-returns backbone: fixes hd and ud ----
  haz_nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    age_int = age_int,
    min_haz = min_haz,
    tol = 1e-14,
    maxiter = 100,
    verbose = FALSE
  )
  
  hd_nr <- haz_nr %>% dplyr::filter(trans == "hd") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
  ud_nr <- haz_nr %>% dplyr::filter(trans == "ud") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
  hu_nr <- haz_nr %>% dplyr::filter(trans == "hu") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
  
  # ---- extract 6 starting values for hu/uh soft-kink ----
  # accepts:
  #   - a named 6-vector: i_hu, s1_hu, s2_hu, i_uh, s1_uh, s2_uh
  #   - an 8/12 full parameter vector, from which hu/uh pieces are extracted
  th0_full <- as_par_vec(theta_start)
  
  get_first_existing <- function(x, nms) {
    hit <- nms[nms %in% names(x)]
    if (length(hit) == 0) return(NA_real_)
    as.numeric(x[[hit[1]]])
  }
  
  th0 <- c(
    i_hu  = get_first_existing(th0_full, c("i_hu")),
    s1_hu = get_first_existing(th0_full, c("s1_hu", "s_hu")),
    s2_hu = get_first_existing(th0_full, c("s2_hu", "s_hu")),
    i_uh  = get_first_existing(th0_full, c("i_uh")),
    s1_uh = get_first_existing(th0_full, c("s1_uh", "s_uh")),
    s2_uh = get_first_existing(th0_full, c("s2_uh", "s_uh"))
  )
  
  if (any(!is.finite(th0))) {
    stop(
      "theta_start must supply hu/uh soft-kink starts: i_hu, s1_hu, s2_hu, i_uh, s1_uh, s2_uh ",
      "(or an 8/12-vector containing the corresponding hu/uh entries).",
      call. = FALSE
    )
  }
  
  # ---- soft-kink expander for just hu and uh ----
  expand_two_softkink <- function(theta6, age, age_int, pivot_age, shape_width, min_haz) {
    theta6 <- as.numeric(theta6)
    names(theta6) <- c("i_hu","s1_hu","s2_hu","i_uh","s1_uh","s2_uh")
    
    age_mid <- age + age_int / 2
    a0 <- min(age)
    x <- age_mid - a0
    
    softplus <- function(z) log1p(exp(z))
    z  <- (age_mid - pivot_age) / shape_width
    z0 <- (min(age_mid) - pivot_age) / shape_width
    kink_term <- shape_width * (softplus(z) - softplus(z0))
    
    log_hu <- theta6["i_hu"] + theta6["s1_hu"] * x +
      (theta6["s2_hu"] - theta6["s1_hu"]) * kink_term
    
    log_uh <- theta6["i_uh"] + theta6["s1_uh"] * x +
      (theta6["s2_uh"] - theta6["s1_uh"]) * kink_term
    
    n <- length(age)
    
    tibble::tibble(
      age = rep(age, times = 2L),
      trans = rep(c("hu", "uh"), each = n),
      hazard = c(
        pmax(exp(log_hu), min_haz),
        pmax(exp(log_uh), min_haz)
      )
    )
  }
  
  # ---- objective: fixed hd/ud, optimize hu/uh over whole age schedule ----
  objective <- function(par6,
                        age,
                        lx,
                        prev_pt,
                        hd_nr,
                        ud_nr,
                        age_int,
                        pivot_age,
                        shape_width,
                        min_haz,
                        w_lx,
                        w_prev) {
    
    haz_two <- expand_two_softkink(
      theta6 = par6,
      age = age,
      age_int = age_int,
      pivot_age = pivot_age,
      shape_width = shape_width,
      min_haz = min_haz
    )
    
    hu_vec <- haz_two %>% dplyr::filter(trans == "hu") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
    uh_vec <- haz_two %>% dplyr::filter(trans == "uh") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
    
    haz_all <- tibble::tibble(
      age = rep(age, times = 4L),
      trans = rep(c("hu", "hd", "ud", "uh"), each = n),
      hazard = c(hu_vec, hd_nr, ud_nr, uh_vec)
    )
    
    P <- haz_to_probs(haz_all, age = age, age_int = age_int)
    lt <- calculate_lt(P, init = c(H = 1, U = 0), age_int = age_int, terminal = "point")
    
    # exact-age targets: use prevalence_point, not interval prevalence
    lx_pred <- lt$lx
    prev_pred <- lt$prevalence_point
    
    loss_lx <- sum((lx_pred - lx)^2, na.rm = TRUE)
    loss_prev <- sum((prev_pred - prev_pt)^2, na.rm = TRUE)
    
    w_lx * loss_lx + w_prev * loss_prev
  }
  
  # ---- simple bounds consistent with your usual sign conventions ----
  lower <- c(
    i_hu  = -Inf,
    s1_hu =  slope_floor,
    s2_hu =  slope_floor,
    i_uh  = -Inf,
    s1_uh = -Inf,
    s2_uh = -Inf
  )
  
  upper <- c(
    i_hu  = Inf,
    s1_hu = Inf,
    s2_hu = Inf,
    i_uh  = Inf,
    s1_uh = -slope_floor,
    s2_uh = -slope_floor
  )
  
  opt <- optim(
    par = th0,
    fn = objective,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = maxit, factr = factr),
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    hd_nr = hd_nr,
    ud_nr = ud_nr,
    age_int = age_int,
    pivot_age = pivot_age,
    shape_width = shape_width,
    min_haz = min_haz,
    w_lx = w_lx,
    w_prev = w_prev
  )
  
  if (verbose) {
    message("derive_returns_hazard_starts(): optim convergence = ", opt$convergence,
            "; objective = ", signif(opt$value, 6))
  }
  
  haz_two <- expand_two_softkink(
    theta6 = opt$par,
    age = age,
    age_int = age_int,
    pivot_age = pivot_age,
    shape_width = shape_width,
    min_haz = min_haz
  )
  
  hu_vec <- haz_two %>% dplyr::filter(trans == "hu") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
  uh_vec <- haz_two %>% dplyr::filter(trans == "uh") %>% dplyr::arrange(age) %>% dplyr::pull(hazard)
  
  haz_out <- tibble::tibble(
    age = rep(age, 4),
    trans = rep(c("hu", "hd", "ud", "uh"), each = n),
    hazard = c(hu_vec, hd_nr, ud_nr, uh_vec)
  )
  
  attr(haz_out, "theta_start_6") <- th0
  attr(haz_out, "theta_hat_6") <- opt$par
  attr(haz_out, "optim") <- opt
  
  haz_out
}
# -----------------------------------------------------------------------------
# These functions replace derive_returns_hazards(),
# but note I'd still like to work in haz_returns_starts()
# creation and Rx df creation in this workflow, i.e.
# that the incoming data should be ground_summary and
# theta_list, (I think); the same tol arg can be used for everything

derive_returns_hazards_precorrection <- function(
    age,
    lx,
    prev_pt,
    Rx,
    hu_target = NULL,
    uh_target = NULL,
    age_int = 1,
    anchor_age = 70,
    min_haz = 1e-12,
    tol = 1e-14,
    maxiter = 100,
    hd_bracket = c(1e-12, 5),
    uh_bracket = c(1e-12, 5),
    diag_tol = 1e-12,
    world_id = NULL,
    verbose = TRUE
) {
  
  age <- as.numeric(age)
  lx <- as.numeric(lx)
  prev_pt <- as.numeric(prev_pt)
  n <- length(age)
  
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  if (length(lx) != n) stop("lx must match age.", call. = FALSE)
  if (length(prev_pt) != n) stop("prev_pt must match age.", call. = FALSE)
  
  Rx <- if (length(Rx) == 1) rep(Rx, n) else as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or n.", call. = FALSE)
  Rx <- pmax(Rx, min_haz)
  
  nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    age_int = age_int,
    min_haz = min_haz,
    tol = tol,
    maxiter = maxiter,
    verbose = FALSE
  )
  
  hu_nr <- nr |> dplyr::filter(trans == "hu") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  hd_nr <- nr |> dplyr::filter(trans == "hd") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  
  hu <- if (is.null(hu_target)) hu_nr else pmax(as.numeric(hu_target), min_haz)
  if (length(hu) != n) stop("hu_target must match age.", call. = FALSE)
  
  if (!is.null(uh_target)) {
    uh_target <- pmax(as.numeric(uh_target), min_haz)
    if (length(uh_target) != n) stop("uh_target must match age.", call. = FALSE)
  }
  
  states <- build_states_targets(age, lx, prev_pt)
  
  hd <- rep(NA_real_, n)
  uh <- rep(NA_real_, n)
  ud <- rep(NA_real_, n)
  
  sanitize_seed <- function(seed, fallback) {
    seed <- as.numeric(seed)[1]
    fallback <- as.numeric(fallback)[1]
    if (!is.finite(seed) || seed <= 0) seed <- fallback
    if (!is.finite(seed) || seed <= 0) seed <- min_haz
    seed
  }
  
  safe_step <- function(H0, U0, hu, hd, uh, ud) {
    if (any(!is.finite(c(H0,U0,hu,hd,uh,ud)))) return(c(H=NA,U=NA))
    Q <- tryCatch(make_Q_returns(hu, hd, uh, ud), error=function(e) NULL)
    if (is.null(Q)) return(c(H=NA,U=NA))
    P <- tryCatch(expm::expm(Q * age_int), error=function(e) NULL)
    if (is.null(P)) return(c(H=NA,U=NA))
    v <- as.numeric(c(H0,U0,0) %*% P)
    c(H=v[1], U=v[2])
  }
  
  solve_hd <- function(hu, uh, Rx, H0, U0, S1, seed) {
    seed <- sanitize_seed(seed, hd_bracket[1])
    
    f <- function(hd) {
      ud <- Rx * hd
      nxt <- safe_step(H0,U0,hu,hd,uh,ud)
      (nxt["H"] + nxt["U"]) - S1
    }
    
    tryCatch(
      uniroot(f, c(min_haz,5), tol=tol, maxiter=maxiter)$root,
      error=function(e) NA_real_
    )
  }
  
  solve_exact <- function(i, hd_seed, uh_seed) {
    hd_seed <- sanitize_seed(hd_seed, hd_nr[i])
    uh_seed <- sanitize_seed(uh_seed, if(!is.null(uh_target)) uh_target[i] else min_haz)
    
    H0 <- states$lH[i]
    U0 <- states$lU[i]
    S1 <- sum(states$nextHU[i,])
    U1 <- states$nextHU[i,2]
    
    f <- function(uh_val) {
      hd_val <- solve_hd(hu[i], uh_val, Rx[i], H0, U0, S1, hd_seed)
      nxt <- safe_step(H0,U0,hu[i],hd_val,uh_val,Rx[i]*hd_val)
      nxt["U"] - U1
    }
    
    uh_val <- tryCatch(
      uniroot(f, c(min_haz,5), tol=tol, maxiter=maxiter)$root,
      error=function(e) NA_real_
    )
    
    if (!is.finite(uh_val)) return(list(hd=NA,uh=NA))
    
    hd_val <- solve_hd(hu[i], uh_val, Rx[i], H0, U0, S1, hd_seed)
    if (!is.finite(hd_val)) return(list(hd=NA,uh=NA))
    
    list(hd=hd_val, uh=uh_val)
  }
  
  make_smooth_uh_ref <- function(i) {
    idx <- seq.int(i + 1, min(i + 3, n - 1))
    idx <- idx[is.finite(uh[idx]) & uh[idx] > 0]
    
    if (length(idx) >= 2) {
      df_fit <- data.frame(
        age = age[idx],
        log_uh = log(uh[idx])
      )
      fit <- stats::lm(log_uh ~ age, data = df_fit)
      ref <- as.numeric(exp(stats::predict(
        fit,
        newdata = data.frame(age = age[i])
      )))[1]
    } else if (i + 1 <= n - 1 && is.finite(uh[i + 1]) && uh[i + 1] > 0) {
      ref <- uh[i + 1]
    } else if (!is.null(uh_target)) {
      ref <- uh_target[i]
    } else {
      ref <- min_haz
    }
    
    pmax(ref, min_haz)
  }
  # --- anchor ---
  anchor <- which.min(abs(age - anchor_age))
  anchor <- max(1, min(anchor, n-1))
  
  sol <- solve_exact(anchor, hd_nr[anchor], uh_target[anchor])
  if (is.finite(sol$hd)) {
    hd[anchor] <- sol$hd
    uh[anchor] <- sol$uh
    ud[anchor] <- Rx[anchor]*hd[anchor]
  } else {
    hd[anchor] <- hd_nr[anchor]
    uh[anchor] <- uh_target[anchor]
    ud[anchor] <- Rx[anchor]*hd[anchor]
  }
  
  # --- forward ---
  for (i in seq(anchor+1, n-1)) {
    sol <- solve_exact(i, hd[i-1], uh[i-1])
    if (is.finite(sol$hd)) {
      hd[i] <- sol$hd
      uh[i] <- sol$uh
    } else {
      hd[i] <- hd[i-1]
      uh[i] <- uh[i-1]
    }
    ud[i] <- Rx[i]*hd[i]
  }
  
  # --- backward ---
  for (i in seq(anchor-1,1,-1)) {
    sol <- solve_exact(i, hd[i+1], uh[i+1])
    if (is.finite(sol$hd)) {
      hd[i] <- sol$hd
      uh[i] <- sol$uh
    } else {
      uh_ref <- make_smooth_uh_ref(i)
      hd_ref <- solve_hd(hu[i], uh_ref, Rx[i],
                         states$lH[i], states$lU[i],
                         sum(states$nextHU[i,]), hd[i+1])
      
      if (is.finite(hd_ref)) {
        hd[i] <- hd_ref
        uh[i] <- uh_ref
      } else {
        hd[i] <- hd[i+1]
        uh[i] <- uh[i+1]
      }
    }
    ud[i] <- Rx[i]*hd[i]
  }
  
  dplyr::bind_rows(
    tibble::tibble(age=age, trans="hu", hazard=hu),
    tibble::tibble(age=age, trans="hd", hazard=hd),
    tibble::tibble(age=age, trans="ud", hazard=ud),
    tibble::tibble(age=age, trans="uh", hazard=uh)
  )
}
# -----------------------------------------------------------------------------


diagnose_returns_branch <- function(
    age,
    lx,
    prev_pt,
    Rx,
    haz,
    age_int = 1,
    age_max = 55,
    hu_bounds_mult = c(0.7, 2.0),
    n_bracket = 21,
    n_refine = 20,
    min_haz = 1e-12,
    tol = 1e-12,
    maxiter = 1000
) {
  
  n <- length(age)
  
  hu0 <- haz |> dplyr::filter(trans=="hu") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  uh0 <- haz |> dplyr::filter(trans=="uh") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  hd0 <- haz |> dplyr::filter(trans=="hd") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  
  states <- build_states_targets(age, lx, prev_pt)
  
  safe_step <- function(H0,U0,hu,hd,uh,ud){
    Q <- make_Q_returns(hu,hd,uh,ud)
    P <- expm::expm(Q * age_int)
    v <- as.numeric(c(H0,U0,0) %*% P)
    c(H=v[1],U=v[2])
  }
  
  solve_hd <- function(hu,uh,Rx,H0,U0,S1,seed){
    f <- function(hd){
      ud <- Rx*hd
      nxt <- safe_step(H0,U0,hu,hd,uh,ud)
      (nxt["H"]+nxt["U"]) - S1
    }
    tryCatch(uniroot(f,c(min_haz,5),tol=tol,maxiter=maxiter)$root,
             error=function(e) NA_real_)
  }
  
  solve_uh <- function(i,hu){
    H0 <- states$lH[i]; U0 <- states$lU[i]
    H1 <- states$nextHU[i,1]; U1 <- states$nextHU[i,2]
    S1 <- H1 + U1
    
    f <- function(uh){
      hd <- solve_hd(hu,uh,Rx[i],H0,U0,S1,hd0[i])
      nxt <- safe_step(H0,U0,hu,hd,uh,Rx[i]*hd)
      nxt["U"] - U1
    }
    
    tryCatch(uniroot(f,c(min_haz,5),tol=tol,maxiter=maxiter)$root,
             error=function(e) NA_real_)
  }
  
  # --- FAST feasible check ---
  feasible <- function(i,hu){
    uh <- solve_uh(i,hu)
    is.finite(uh)
  }
  
  ages_use <- which(age <= age_max & seq_len(n) < n)
  
  res <- lapply(ages_use, function(i){
    
    lo <- hu0[i] * hu_bounds_mult[1]
    hi <- hu0[i] * hu_bounds_mult[2]
    
    # --- coarse bracket scan ---
    grid <- seq(lo, hi, length.out = n_bracket)
    ok <- sapply(grid, function(h) feasible(i,h))
    
    if (!any(ok)) {
      return(data.frame(
        age = age[i],
        hu_min_feasible = NA,
        hu_max_feasible = NA
      ))
    }
    
    # first feasible index
    j <- which(ok)[1]
    
    # bracket
    lo2 <- if (j == 1) grid[1] else grid[j-1]
    hi2 <- grid[j]
    
    # --- refine (bisection-style) ---
    for (k in seq_len(n_refine)) {
      mid <- 0.5 * (lo2 + hi2)
      if (feasible(i, mid)) {
        hi2 <- mid
      } else {
        lo2 <- mid
      }
    }
    
    data.frame(
      age = age[i],
      hu_min_feasible = hi2,
      hu_max_feasible = max(grid[ok])
    )
  })
  
  dplyr::bind_rows(res)
}



correct_returns_hazards <- function(
    pre_haz,
    diag_info,
    age,
    lx,
    prev_pt,
    Rx,
    age_int = 1,
    correction_max_age = 55,
    min_haz = 1e-12,
    tol = 1e-14,
    maxiter = 200,
    eps = 0.01,
    verbose = TRUE
) {
  
  n <- length(age)
  
  hu_pre <- pre_haz |>
    dplyr::filter(trans == "hu") |>
    dplyr::arrange(age) |>
    dplyr::pull(hazard)
  
  hd_pre <- pre_haz |>
    dplyr::filter(trans == "hd") |>
    dplyr::arrange(age) |>
    dplyr::pull(hazard)
  
  uh_pre <- pre_haz |>
    dplyr::filter(trans == "uh") |>
    dplyr::arrange(age) |>
    dplyr::pull(hazard)
  
  # full-length floor vector; outside correction ages, leave as hu_pre
  hu_floor_full <- hu_pre
  idx_diag <- match(diag_info$age, age)
  idx_diag <- idx_diag[is.finite(idx_diag)]
  hu_floor_full[idx_diag] <- diag_info$hu_min_feasible[match(age[idx_diag], diag_info$age)]
  
  # correction window
  idx_corr <- which(age <= correction_max_age)
  if (length(idx_corr) == 0) {
    return(pre_haz)
  }
  
  i_left <- min(idx_corr)
  i_right <- max(idx_corr)
  
  hu_left  <- hu_floor_full[i_left] * (1 + eps)
  hu_right <- hu_pre[i_right]
  
  hu_new <- hu_pre
  
  if (length(idx_corr) == 1) {
    hu_new[idx_corr] <- hu_left
  } else {
    x <- age[idx_corr]
    log_left <- log(hu_left)
    log_right <- log(hu_right)
    log_bridge <- log_left + (log_right - log_left) *
      (x - x[1]) / (x[length(x)] - x[1])
    hu_new[idx_corr] <- exp(log_bridge)
  }
  
  # enforce feasibility ONLY on correction ages
  hu_new[idx_corr] <- pmax(
    hu_new[idx_corr],
    hu_floor_full[idx_corr] * (1 + 1e-8)
  )
  hu_new <- pmax(hu_new, min_haz)
  
  states <- build_states_targets(age, lx, prev_pt)
  
  safe_step <- function(H0, U0, hu, hd, uh, ud) {
    vals <- c(H0, U0, hu, hd, uh, ud)
    if (any(!is.finite(vals)) || any(vals[c(3,4,5,6)] < 0)) {
      return(c(H = NA_real_, U = NA_real_))
    }
    Q <- tryCatch(make_Q_returns(hu, hd, uh, ud), error = function(e) NULL)
    if (is.null(Q)) return(c(H = NA_real_, U = NA_real_))
    P <- tryCatch(expm::expm(Q * age_int), error = function(e) NULL)
    if (is.null(P) || any(!is.finite(P))) return(c(H = NA_real_, U = NA_real_))
    v <- as.numeric(c(H0, U0, 0) %*% P)
    c(H = v[1], U = v[2])
  }
  
  solve_hd <- function(hu, uh, Rx_i, H0, U0, S1, seed) {
    seed <- if (is.finite(seed) && seed > 0) seed else min_haz
    
    f <- function(hd) {
      if (!is.finite(hd) || hd <= 0) return(NA_real_)
      ud <- Rx_i * hd
      nxt <- safe_step(H0, U0, hu, hd, uh, ud)
      if (any(!is.finite(nxt))) return(NA_real_)
      (nxt["H"] + nxt["U"]) - S1
    }
    
    lo <- max(min_haz, seed / 20)
    hi <- max(10 * seed, 1)
    
    f_lo <- f(lo)
    f_hi <- f(hi)
    
    k <- 0L
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 40L) {
      hi <- hi * 2
      f_hi <- f(hi)
      k <- k + 1L
    }
    
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      return(NA_real_)
    }
    
    tryCatch(
      uniroot(f, c(lo, hi), tol = tol, maxiter = maxiter)$root,
      error = function(e) NA_real_
    )
  }
  
  solve_exact <- function(i, hd_seed, uh_seed) {
    H0 <- states$lH[i]
    U0 <- states$lU[i]
    S1 <- sum(states$nextHU[i, ])
    U1 <- states$nextHU[i, 2]
    
    hd_seed <- if (is.finite(hd_seed) && hd_seed > 0) hd_seed else hd_pre[i]
    uh_seed <- if (is.finite(uh_seed) && uh_seed > 0) uh_seed else uh_pre[i]
    if (!is.finite(uh_seed) || uh_seed <= 0) uh_seed <- min_haz
    
    f <- function(uh_val) {
      if (!is.finite(uh_val) || uh_val <= 0) return(NA_real_)
      hd_val <- solve_hd(hu_new[i], uh_val, Rx[i], H0, U0, S1, hd_seed)
      if (!is.finite(hd_val)) return(NA_real_)
      nxt <- safe_step(H0, U0, hu_new[i], hd_val, uh_val, Rx[i] * hd_val)
      if (any(!is.finite(nxt))) return(NA_real_)
      nxt["U"] - U1
    }
    
    lo <- max(min_haz, uh_seed / 50)
    hi <- max(50 * uh_seed, 1)
    
    f_lo <- f(lo)
    f_hi <- f(hi)
    
    k <- 0L
    while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 40L) {
      hi <- hi * 2
      f_hi <- f(hi)
      k <- k + 1L
    }
    
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      return(list(ok = FALSE))
    }
    
    uh_val <- tryCatch(
      uniroot(f, c(lo, hi), tol = tol, maxiter = maxiter)$root,
      error = function(e) NA_real_
    )
    if (!is.finite(uh_val)) return(list(ok = FALSE))
    
    hd_val <- solve_hd(hu_new[i], uh_val, Rx[i], H0, U0, S1, hd_seed)
    if (!is.finite(hd_val)) return(list(ok = FALSE))
    
    list(ok = TRUE, hd = hd_val, uh = uh_val)
  }
  
  hd <- uh <- ud <- rep(NA_real_, n)
  
  anchor <- which.min(abs(age - 70))
  anchor <- max(1, min(anchor, n - 1))
  
  sol <- solve_exact(anchor, hd_pre[anchor], uh_pre[anchor])
  if (!sol$ok) stop("Anchor solve failed in correction.")
  
  hd[anchor] <- sol$hd
  uh[anchor] <- sol$uh
  ud[anchor] <- Rx[anchor] * hd[anchor]
  
  for (i in seq.int(anchor + 1, n - 1)) {
    sol <- solve_exact(i, hd[i - 1], uh[i - 1])
    if (!sol$ok) stop(sprintf("Forward solve failed at age %s", age[i]))
    hd[i] <- sol$hd
    uh[i] <- sol$uh
    ud[i] <- Rx[i] * hd[i]
  }
  
  for (i in seq.int(anchor - 1, 1, by = -1)) {
    sol <- solve_exact(i, hd[i + 1], uh[i + 1])
    if (!sol$ok) stop(sprintf("Backward solve failed at age %s", age[i]))
    hd[i] <- sol$hd
    uh[i] <- sol$uh
    ud[i] <- Rx[i] * hd[i]
  }
  
  extrap <- function(x) {
    if (n >= 3 && all(is.finite(x[(n - 2):(n - 1)])) && all(x[(n - 2):(n - 1)] > 0)) {
      exp(2 * log(x[n - 1]) - log(x[n - 2]))
    } else {
      x[n - 1]
    }
  }
  
  hd[n] <- extrap(hd)
  uh[n] <- extrap(uh)
  ud[n] <- Rx[n] * hd[n]
  
  dplyr::bind_rows(
    tibble::tibble(age = age, trans = "hu", hazard = hu_new),
    tibble::tibble(age = age, trans = "hd", hazard = hd),
    tibble::tibble(age = age, trans = "ud", hazard = ud),
    tibble::tibble(age = age, trans = "uh", hazard = uh)
  )
}


derive_returns_hazards <- function(
    age,
    lx,
    prev_pt,
    Rx,
    theta_start,
    age_int = 1,
    pivot_age = 75,
    shape_width = 5,
    anchor_age = 70,
    correction_max_age = 55,
    min_haz = 1e-12,
    tol = 1e-14,
    maxiter = 200,
    factr = 1e6,
    w_lx = 1,
    w_prev = 1,
    slope_floor = 1e-6,
    hd_bracket = c(1e-12, 5),
    uh_bracket = c(1e-12, 5),
    diag_tol = 1e-12,
    hu_bounds_mult = c(0.7, 2.0),
    n_bracket = 21,
    n_refine = 20,
    eps = 0.01,
    world_id = NULL,
    verbose = TRUE
) {
  
  age <- as.numeric(age)
  lx <- as.numeric(lx)
  prev_pt <- as.numeric(prev_pt)
  n <- length(age)
  
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  if (length(lx) != n) stop("lx must match age.", call. = FALSE)
  if (length(prev_pt) != n) stop("prev_pt must match age.", call. = FALSE)
  
  if (length(Rx) == 1) Rx <- rep(as.numeric(Rx), n)
  Rx <- as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or n.", call. = FALSE)
  Rx <- pmax(Rx, min_haz)
  
  # ---------------------------------------------------------------------------
  # Step 1: derive internal starts for hu / uh
  # ---------------------------------------------------------------------------
  hz_start <- derive_returns_hazard_starts(
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    theta_start = theta_start,
    age_int = age_int,
    pivot_age = pivot_age,
    shape_width = shape_width,
    min_haz = min_haz,
    maxit = maxiter,
    factr = factr,
    w_lx = w_lx,
    w_prev = w_prev,
    slope_floor = slope_floor,
    verbose = verbose
  )
  
  hu_target <- hz_start |>
    dplyr::filter(.data$trans == "hu") |>
    dplyr::arrange(.data$age) |>
    dplyr::pull(.data$hazard)
  
  uh_target <- hz_start |>
    dplyr::filter(.data$trans == "uh") |>
    dplyr::arrange(.data$age) |>
    dplyr::pull(.data$hazard)
  
  if (length(hu_target) != n) stop("Internal hu_target length mismatch.", call. = FALSE)
  if (length(uh_target) != n) stop("Internal uh_target length mismatch.", call. = FALSE)
  
  # ---------------------------------------------------------------------------
  # Step 2: precorrection
  # ---------------------------------------------------------------------------
  pre <- derive_returns_hazards_precorrection(
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    hu_target = hu_target,
    uh_target = uh_target,
    age_int = age_int,
    anchor_age = anchor_age,
    min_haz = min_haz,
    tol = tol,
    maxiter = maxiter,
    hd_bracket = hd_bracket,
    uh_bracket = uh_bracket,
    diag_tol = diag_tol,
    world_id = world_id,
    verbose = verbose
  )
  
  # ---------------------------------------------------------------------------
  # Step 3: diagnose early-age feasible branch
  # ---------------------------------------------------------------------------
  diagi <- diagnose_returns_branch(
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    haz = pre,
    age_int = age_int,
    age_max = correction_max_age,
    hu_bounds_mult = hu_bounds_mult,
    n_bracket = n_bracket,
    n_refine = n_refine,
    min_haz = min_haz,
    tol = diag_tol,
    maxiter = maxiter
  )
  
  # ---------------------------------------------------------------------------
  # Step 4: correct
  # ---------------------------------------------------------------------------
  final <- correct_returns_hazards(
    pre_haz = pre,
    diag_info = diagi,
    age = age,
    lx = lx,
    prev_pt = prev_pt,
    Rx = Rx,
    age_int = age_int,
    correction_max_age = correction_max_age,
    min_haz = min_haz,
    tol = tol,
    maxiter = maxiter,
    eps = eps,
    verbose = verbose
  )
  
  final
}


haz_block_to_lt <- function(hz_block, init, age_int) {
  P <- haz_to_probs(
    hazard = hz_block %>% select(age, trans, hazard),
    age = age,
    age_int = age_int
  )
  
  calculate_lt(P, init = init, age_int = age_int)
}
