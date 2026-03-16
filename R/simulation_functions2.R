# simulation_functions.R
# Clean, harmonized function set for decomp_sullivan multistate simulation workflow.
# - Canonical make_Q() used by haz_to_probs() is preserved.
# - Candidate world generation via SVD near-null directions.
# - Single-pass polishing via polish_many() (no second-stage refinement).
# - Exact "returns" hazards obtained by snapping (hu0,uh0) to per-age isoclines.
# - Deterministic no-returns hazards via derive_noreturns_hazards().
#
# Dependencies: dplyr, tidyr, purrr, tibble, expm

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(expm)
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
# calculate_lt <- function(P, init = c(H = 1, U = 0), age_int = 1) {
#   P |>
#     tidyr::pivot_longer(-c(age, from), names_to = "to", values_to = "p") |>
#     dplyr::filter(.data$from != "D") |>
#     tidyr::unite("state", c("from","to"), sep = "") |>
#     tidyr::pivot_wider(names_from = state, values_from = p) |>
#     dplyr::group_modify(~{
#       .x$lh <- numeric(nrow(.x))
#       .x$lu <- numeric(nrow(.x))
#       
#       init_use <- init
#       .x$lh[1] <- init_use["H"]
#       .x$lu[1] <- init_use["U"]
#       
#       for (ii in seq_len(nrow(.x) - 1)) {
#         .x$lh[ii + 1] <- .x$lh[ii] * .x$HH[ii] + .x$lu[ii] * .x$UH[ii]
#         .x$lu[ii + 1] <- .x$lu[ii] * .x$UU[ii] + .x$lh[ii] * .x$HU[ii]
#       }
#       .x
#     }) |>
#     dplyr::mutate(
#       prevalence = .data$lu / (.data$lh + .data$lu),
#       lx         = .data$lh + .data$lu,
#       Lx         = age_int * (.data$lx + dplyr::lead(.data$lx, default = 0)) / 2
#     )
# }
# Experimental override for derive_noreturns_hazards()
# Uses interval-average unhealthy share (trapezoid proxy) instead of point prevalence
# for constructing the target unhealthy stock vector inside the no-returns solver.
#
# Source AFTER simulation_functions.R


# ------------------------------------------------------------------
# Prevalence handling
# ------------------------------------------------------------------
# calculate_lt() now returns:
#   prevalence_point    = lu / lx at exact age x
#   prevalence_interval = Lu / Lx over the interval [x, x+age_int)
# and sets prevalence = prevalence_interval by default.
#
# Consistency rule used below:
#   - no-returns calibration works off the prevalence vector supplied by caller;
#     by default that is treated as already being on the intended scale ("as_is").
#   - state-target reconstruction for returns matching uses point prevalence,
#     because lh_x and lu_x are exact-age stocks.
# ------------------------------------------------------------------

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

build_states_targets <- function(age, lx, prev, prevalence_point = NULL) {
  age <- as.numeric(age)
  lx  <- as.numeric(lx)
  n   <- length(age)
  
  if (length(lx) != n) {
    stop("lx must have same length as age.", call. = FALSE)
  }
  
  prev_state <- .coerce_state_prev(
    prev = prev,
    prevalence_point = prevalence_point,
    lx = lx,
    n = n
  )
  prev_state <- pmin(pmax(prev_state, 0), 1)
  
  lH <- lx * (1 - prev_state)
  lU <- lx * prev_state
  lD <- 1 - lx
  
  nextHU <- cbind(lH, lU)
  if (n >= 2) {
    nextHU[1:(n - 1), ] <- cbind(lH[2:n], lU[2:n])
    nextHU[n, ] <- nextHU[n - 1, ]
  }
  
  list(
    age = age,
    lx = lx,
    prev = prev_state,
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

derive_noreturns_hazards <- function(age,
                                     lx,
                                     prev,
                                     Rx,
                                     age_int = 1,
                                     min_haz = 1e-12,
                                     tol = 1e-20,
                                     maxiter = 100,
                                     hu_bracket = c(1e-12, 5),
                                     hd_bracket = c(1e-12, 5),
                                     verbose = TRUE,
                                     prev_mode = c("as_is", "point_to_interval"),
                                     interval_last = c("carry", "point")) {
  
  prev_mode <- match.arg(prev_mode)
  interval_last <- match.arg(interval_last)
  
  age <- as.numeric(age)
  n <- length(age)
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  
  lx <- as.numeric(lx)
  prev <- as.numeric(prev)
  if (length(lx) != n) stop("lx must have same length as age.", call. = FALSE)
  if (length(prev) != n) stop("prev must have same length as age.", call. = FALSE)
  
  if (length(Rx) == 1) Rx <- rep(as.numeric(Rx), n)
  Rx <- as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or length(age).", call. = FALSE)
  
  prev_used <- if (prev_mode == "point_to_interval") {
    interval_prev_trapezoid(lx = lx, prev = prev, last = interval_last)
  } else {
    prev
  }
  
  # Look here!
  lu_tgt <- prev_used * lx
  lh_tgt <- lx - lu_tgt
  
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

derive_returns_hazards_from_Rx <- function(
    age,
    lx,
    prev,
    Rx,
    age_int = 1,
    bounds = c(1e-12, 2),
    hu_bounds = NULL,
    uh_bounds = NULL,
    maxit = 500,
    reltol = 1e-14,
    hu_0,
    uh_0,
    smooth_w = 0.5,
    curvature_w = 0.5,
    bound_w  = 1e-3,
    turnover_K = 1.2,
    turnover_w = 1e4,
    fit_w = 1e6,
    solver = c("project"),
    snap_tol = 1e-20,
    verbose = FALSE,
    eps_log = 1e-2,
    tiny = 1e-14,
    line_tol = 1e-12,
    final_refine = TRUE,
    refine_tol = 1e-16,
    refine_dist_w = 1e-6,
    extrap_last_age = TRUE,
    prevalence_point = NULL
){
  solver <- match.arg(solver)
  
  age <- as.numeric(age)
  lx <- as.numeric(lx)
  prev <- as.numeric(prev)
  Rx <- as.numeric(Rx)
  hu_0 <- as.numeric(hu_0)
  uh_0 <- as.numeric(uh_0)
  
  n <- length(age)
  stopifnot(length(lx) == n, length(prev) == n, length(Rx) == n,
            length(hu_0) == n, length(uh_0) == n)
  
  hz_nr <- derive_noreturns_hazards(
    age = age,
    lx = lx,
    prev = prev,
    Rx = Rx,
    age_int = age_int,
    verbose = FALSE
  )
  
  hd <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
  ud <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
  hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
  
  if (is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  if (is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
  
  hu_lo <- pmax(hu_bounds[, 1], hu_nr)
  hu_hi <- hu_bounds[, 2]
  uh_lo <- pmax(0, uh_bounds[, 1])
  uh_hi <- uh_bounds[, 2]
  
  states <- build_states_targets(age, lx, prev, prevalence_point = prevalence_point)
  
  hu_out <- numeric(n)
  uh_out <- numeric(n)
  p2_loss <- rep(NA_real_, n)
  degenerate <- rep(FALSE, n)
  
  for (i in seq_len(n)) {
    cap <- turnover_K * (hu_0[i] + uh_0[i])
    
    p1 <- c(clamp(hu_nr[i], hu_lo[i], hu_hi[i]), 0)
    p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
    
    p0 <- c(clamp(hu_0[i], hu_lo[i], hu_hi[i]), clamp(uh_0[i], uh_lo[i], uh_hi[i]))
    
    hu2_try <- clamp(p1[1] * exp(eps_log), hu_lo[i], hu_hi[i])
    uh_hi_eff <- min(uh_hi[i], cap - hu2_try)
    if (!is.finite(uh_hi_eff)) uh_hi_eff <- uh_hi[i]
    uh_hi_eff <- max(uh_hi_eff, uh_lo[i])
    
    sol2 <- solve_uh_given_hu(
      i, states, hd, ud,
      hu = hu2_try,
      uh_lo = uh_lo[i], uh_hi = uh_hi_eff,
      age_int = age_int, tiny = tiny
    )
    
    have_p2 <- is.finite(sol2$loss) && is.finite(sol2$uh) &&
      (abs(hu2_try - p1[1]) > 0) && (sol2$loss <= line_tol)
    
    if (!have_p2) {
      uh_seed <- clamp(max(uh_lo[i], tiny) * exp(eps_log), uh_lo[i], uh_hi[i])
      hu_hi_eff <- min(hu_hi[i], cap - uh_seed)
      if (!is.finite(hu_hi_eff)) hu_hi_eff <- hu_hi[i]
      hu_hi_eff <- max(hu_hi_eff, hu_lo[i])
      
      sol2b <- solve_hu_given_uh(
        i, states, hd, ud,
        uh = uh_seed,
        hu_lo = hu_lo[i], hu_hi = hu_hi_eff,
        age_int = age_int
      )
      
      have_p2 <- is.finite(sol2b$loss) && is.finite(sol2b$hu) &&
        (abs(uh_seed - p1[2]) > 0) && (sol2b$loss <= line_tol)
      
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
        
        if (isTRUE(final_refine)) {
          f_t <- function(t) {
            pt <- p1 + t * v
            ls <- loss_one_age(i, states, hd, ud, hu = pt[1], uh = pt[2], age_int = age_int)
            ls + refine_dist_w * (t - t_proj)^2
          }
          opt <- optimize(f_t, interval = c(feas$t_lo, feas$t_hi), tol = refine_tol)
          t_star <- opt$minimum
        } else {
          t_star <- t_proj
        }
        
        p_star <- p1 + t_star * v
        hu_hat <- p_star[1]
        uh_hat <- p_star[2]
      }
    }
    
    hu_out[i] <- hu_hat
    uh_out[i] <- uh_hat
    
    if (verbose && (i %% 10 == 0 || i == 1 || i == n)) {
      ls <- loss_one_age(i, states, hd, ud, hu = hu_out[i], uh = uh_out[i], age_int = age_int)
      message(sprintf("age %s: loss=%.3e  p2_loss=%s  deg=%s",
                      age[i], ls, format(p2_loss[i], scientific = TRUE, digits = 2), degenerate[i]))
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
  
  attr(out, "diagnostics") <- list(p2_loss = p2_loss, degenerate = degenerate)
  out
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

# ---- Exact returns snapper (isocline projection) ----

resid_one_age <- function(i, states, hd, ud, hu, uh, age_int = 1){
  cur <- c(states$lH[i], states$lU[i], states$lD[i])
  targ <- c(states$nextHU[i,1], states$nextHU[i,2])
  Q <- make_Q_returns(hu=hu, hd=hd[i], uh=uh, ud=ud[i])
  P <- expm::expm(Q * age_int)
  nxt <- as.numeric(cur %*% P)
  c(rH = nxt[1] - targ[1], rU = nxt[2] - targ[2])
}

loss_one_age <- function(i, states, hd, ud, hu, uh, age_int = 1){
  r <- resid_one_age(i, states, hd, ud, hu, uh, age_int)
  sum(r*r)
}

solve_uh_given_hu <- function(i, states, hd, ud, hu,
                              uh_lo, uh_hi,
                              age_int = 1,
                              tiny = 1e-14){
  if(!is.finite(hu) || hu <= 0) return(list(uh=NA_real_, loss=Inf))
  if(uh_hi < uh_lo) return(list(uh=NA_real_, loss=Inf))
  
  f <- function(z){
    uh <- exp(z) - tiny
    uh <- clamp(uh, uh_lo, uh_hi)
    loss_one_age(i, states, hd, ud, hu=hu, uh=uh, age_int=age_int)
  }
  zlo <- log(uh_lo + tiny)
  zhi <- log(uh_hi + tiny)
  opt <- optimize(f, interval = c(zlo, zhi), tol = 1e-20)
  uh_hat <- clamp(exp(opt$minimum) - tiny, uh_lo, uh_hi)
  list(uh = uh_hat, loss = opt$objective)
}

solve_hu_given_uh <- function(i, states, hd, ud, uh,
                              hu_lo, hu_hi,
                              age_int = 1){
  if(!is.finite(uh) || uh < 0) return(list(hu=NA_real_, loss=Inf))
  if(hu_hi < hu_lo) return(list(hu=NA_real_, loss=Inf))
  
  f <- function(loghu){
    hu <- exp(loghu)
    loss_one_age(i, states, hd, ud, hu=hu, uh=uh, age_int=age_int)
  }
  lo <- log(pmax(hu_lo, 1e-300))
  hi <- log(pmax(hu_hi, 1e-300))
  opt <- optimize(f, interval = c(lo, hi), tol = 1e-12)
  hu_hat <- clamp(exp(opt$minimum), hu_lo, hu_hi)
  list(hu = hu_hat, loss = opt$objective)
}

project_point_to_line <- function(p, a, v){
  vv <- sum(v*v)
  if(vv <= 0) return(list(t = 0, proj = a))
  t <- sum((p - a) * v) / vv
  list(t = t, proj = a + t * v)
}

clamp_t_feasible <- function(a, v,
                             hu_lo, hu_hi,
                             uh_lo, uh_hi,
                             cap = Inf){
  t_lo <- -Inf
  t_hi <-  Inf
  
  add_ge <- function(beta, rhs){
    if(abs(beta) < 1e-30){
      if(rhs > 0) return(FALSE) else return(TRUE)
    }
    t0 <- rhs / beta
    if(beta > 0) t_lo <<- max(t_lo, t0) else t_hi <<- min(t_hi, t0)
    TRUE
  }
  add_le <- function(beta, rhs){
    if(abs(beta) < 1e-30){
      if(rhs < 0) return(FALSE) else return(TRUE)
    }
    t0 <- rhs / beta
    if(beta > 0) t_hi <<- min(t_hi, t0) else t_lo <<- max(t_lo, t0)
    TRUE
  }
  
  if(!add_ge(v[1], hu_lo - a[1])) return(list(ok=FALSE))
  if(!add_le(v[1], hu_hi - a[1])) return(list(ok=FALSE))
  if(!add_ge(v[2], uh_lo - a[2])) return(list(ok=FALSE))
  if(!add_le(v[2], uh_hi - a[2])) return(list(ok=FALSE))
  
  if(is.finite(cap)){
    if(!add_le(v[1] + v[2], cap - (a[1] + a[2]))) return(list(ok=FALSE))
  }
  
  ok <- is.finite(t_lo) && is.finite(t_hi) && (t_lo <= t_hi)
  list(ok=ok, t_lo=t_lo, t_hi=t_hi)
}


# derive_returns_hazards_from_Rx <- function(
    #     age,
#     lx,
#     prev,
#     Rx,
#     age_int = 1,
#     bounds = c(1e-12, 2),
#     hu_bounds = NULL,
#     uh_bounds = NULL,
#     maxit = 500,           # accepted for back-compat; unused
#     reltol = 1e-14,        # accepted for back-compat; unused
#     hu_0,
#     uh_0,
#     smooth_w = 0.5,        # accepted; unused
#     curvature_w = 0.5,     # accepted; unused
#     bound_w  = 1e-3,       # accepted; unused
#     turnover_K = 1.2,
#     turnover_w = 1e4,      # accepted; unused
#     fit_w = 1e6,           # accepted; unused
#     solver = c("project"), # accepted; unused (always projection+refine)
#     snap_tol = 1e-20,      # diagnostic only
#     verbose = FALSE,
#     # snapper controls
#     eps_log = 1e-2,
#     tiny = 1e-14,
#     line_tol = 1e-12,
#     # NEW: numerical tighten + last-age behaviour
#     final_refine = TRUE,
#     refine_tol = 1e-16,
#     refine_dist_w = 1e-6,
#     extrap_last_age = TRUE
# ){
#   
#   solver <- match.arg(solver)
#   
#   n <- length(age)
#   stopifnot(length(lx) == n, length(prev) == n, length(Rx) == n,
#             length(hu_0) == n, length(uh_0) == n)
#   
#   # Deterministic no-returns gives hd, ud, and HU lower bound (hu_nr)
#   hz_nr <- derive_noreturns_hazards(
#     age = age,
#     lx = lx,
#     prev = prev,
#     Rx = Rx,
#     age_int = age_int,
#     verbose = FALSE
#   )
#   
#   hd <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
#   ud <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
#   hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
#   
#   # bounds
#   if(is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   if(is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   
#   hu_lo <- pmax(hu_bounds[,1], hu_nr)
#   hu_hi <- hu_bounds[,2]
#   uh_lo <- pmax(0, uh_bounds[,1])
#   uh_hi <- uh_bounds[,2]
#   
#   states <- build_states_targets(age, lx, prev)
#   
#   hu_out <- numeric(n)
#   uh_out <- numeric(n)
#   p2_loss <- rep(NA_real_, n)
#   degenerate <- rep(FALSE, n)
#   
#   for(i in seq_len(n)){
#     
#     cap <- turnover_K * (hu_0[i] + uh_0[i])
#     
#     # P1: exact no-returns point on the isocline
#     p1 <- c(clamp(hu_nr[i], hu_lo[i], hu_hi[i]), 0)
#     p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
#     
#     # P0: stage-1 point to project
#     p0 <- c(clamp(hu_0[i], hu_lo[i], hu_hi[i]), clamp(uh_0[i], uh_lo[i], uh_hi[i]))
#     
#     # Build P2 to define isocline direction (best-effort)
#     hu2_try <- clamp(p1[1] * exp(eps_log), hu_lo[i], hu_hi[i])
#     uh_hi_eff <- min(uh_hi[i], cap - hu2_try)
#     if(!is.finite(uh_hi_eff)) uh_hi_eff <- uh_hi[i]
#     uh_hi_eff <- max(uh_hi_eff, uh_lo[i])
#     
#     sol2 <- solve_uh_given_hu(
#       i, states, hd, ud,
#       hu = hu2_try,
#       uh_lo = uh_lo[i], uh_hi = uh_hi_eff,
#       age_int = age_int, tiny = tiny
#     )
#     
#     have_p2 <- is.finite(sol2$loss) && is.finite(sol2$uh) &&
#       (abs(hu2_try - p1[1]) > 0) && (sol2$loss <= line_tol)
#     
#     if(!have_p2){
#       uh_seed <- clamp(max(uh_lo[i], tiny) * exp(eps_log), uh_lo[i], uh_hi[i])
#       hu_hi_eff <- min(hu_hi[i], cap - uh_seed)
#       if(!is.finite(hu_hi_eff)) hu_hi_eff <- hu_hi[i]
#       hu_hi_eff <- max(hu_hi_eff, hu_lo[i])
#       
#       sol2b <- solve_hu_given_uh(
#         i, states, hd, ud,
#         uh = uh_seed,
#         hu_lo = hu_lo[i], hu_hi = hu_hi_eff,
#         age_int = age_int
#       )
#       
#       have_p2 <- is.finite(sol2b$loss) && is.finite(sol2b$hu) &&
#         (abs(uh_seed - p1[2]) > 0) && (sol2b$loss <= line_tol)
#       
#       if(have_p2){
#         p2 <- c(sol2b$hu, uh_seed)
#         p2_loss[i] <- sol2b$loss
#       } else {
#         # accept best-effort direction if it differs from p1 (avoid degeneracy)
#         if(is.finite(sol2$loss) && is.finite(sol2$uh) && abs(hu2_try - p1[1]) > 0){
#           p2 <- c(hu2_try, sol2$uh)
#           p2_loss[i] <- sol2$loss
#           have_p2 <- TRUE
#         } else if(is.finite(sol2b$loss) && is.finite(sol2b$hu) && abs(uh_seed - p1[2]) > 0){
#           p2 <- c(sol2b$hu, uh_seed)
#           p2_loss[i] <- sol2b$loss
#           have_p2 <- TRUE
#         } else {
#           p2 <- p1
#           have_p2 <- FALSE
#         }
#       }
#     } else {
#       p2 <- c(hu2_try, sol2$uh)
#       p2_loss[i] <- sol2$loss
#     }
#     
#     v <- p2 - p1
#     if(sum(v*v) < 1e-30){
#       degenerate[i] <- TRUE
#       hu_hat <- p1[1]
#       uh_hat <- p1[2]
#     } else {
#       
#       # feasible t interval on the segment (bounds + turnover cap)
#       feas <- clamp_t_feasible(
#         a = p1, v = v,
#         hu_lo = hu_lo[i], hu_hi = hu_hi[i],
#         uh_lo = uh_lo[i], uh_hi = uh_hi[i],
#         cap = cap
#       )
#       
#       if(!isTRUE(feas$ok)){
#         hu_hat <- p1[1]
#         uh_hat <- p1[2]
#       } else {
#         # orthogonal projection parameter t
#         pr <- project_point_to_line(p0, p1, v)
#         t_proj <- clamp(pr$t, feas$t_lo, feas$t_hi)
#         
#         if(isTRUE(final_refine)){
#           # refine along the line segment to remove expm/optimize numerical noise
#           f_t <- function(t){
#             pt <- p1 + t * v
#             ls <- loss_one_age(i, states, hd, ud, hu = pt[1], uh = pt[2], age_int = age_int)
#             ls + refine_dist_w * (t - t_proj)^2
#           }
#           opt <- optimize(f_t, interval = c(feas$t_lo, feas$t_hi), tol = refine_tol)
#           t_star <- opt$minimum
#         } else {
#           t_star <- t_proj
#         }
#         
#         p_star <- p1 + t_star * v
#         hu_hat <- p_star[1]
#         uh_hat <- p_star[2]
#         uh_hat <- solve_uh_given_hu
#       }
#     }
#     
#     hu_out[i] <- hu_hat
#     uh_out[i] <- uh_hat
#     
#     if(verbose && (i %% 10 == 0 || i == 1 || i == n)){
#       ls <- loss_one_age(i, states, hd, ud, hu = hu_out[i], uh = uh_out[i], age_int = age_int)
#       message(sprintf("age %s: loss=%.3e  p2_loss=%s  deg=%s",
#                       age[i], ls, format(p2_loss[i], scientific=TRUE, digits=2), degenerate[i]))
#     }
#   }
#   
#   # Optional: last-age extrapolation for returns hazards too (doesn't affect matching)
#   if(isTRUE(extrap_last_age) && n >= 3){
#     min_haz <- bounds[1]
#     extrap_last <- function(x) {
#       if (is.finite(x[n-1]) && is.finite(x[n-2]) && x[n-1] > 0 && x[n-2] > 0) {
#         exp(2*log(x[n-1]) - log(x[n-2]))
#       } else if (is.finite(x[n-1])) {
#         x[n-1]
#       } else {
#         min_haz
#       }
#     }
#     
#     # Extrapolate and then enforce bounds + turnover cap
#     hu_out[n] <- clamp(pmax(extrap_last(hu_out), min_haz), hu_lo[n], hu_hi[n])
#     uh_out[n] <- clamp(pmax(extrap_last(uh_out), min_haz), uh_lo[n], uh_hi[n])
#     
#     cap_n <- turnover_K * (hu_0[n] + uh_0[n])
#     if(is.finite(cap_n) && (hu_out[n] + uh_out[n] > cap_n)){
#       # keep hu fixed, trim uh
#       uh_out[n] <- max(uh_lo[n], cap_n - hu_out[n])
#     }
#   }
#   
#   out <- tibble::tibble(
#     age = rep(age, times = 4),
#     trans = rep(c("hd","hu","ud","uh"), each = n),
#     hazard = c(hd, hu_out, ud, uh_out)
#   )
#   
#   attr(out, "diagnostics") <- list(p2_loss = p2_loss, degenerate = degenerate)
#   out
# }
# Drop-in replacement for derive_returns_hazards_from_Rx()
#
# Intended use:
#   source("R/simulation_functions.R")
#   source("derive_returns_hazards_from_Rx_inplace_refine.R")
#
# This redefines derive_returns_hazards_from_Rx() using the existing machinery
# already present in simulation_functions.R, with one minimal extra refinement:
# immediately after the line-snap yields (hu_hat, uh_hat), hold hu_hat fixed and
# call solve_uh_given_hu() one more time to put the point on the implied one-age
# curve for that hu.

# Experimental curve-snap alternative to derive_returns_hazards_from_Rx()
#
# Assumes simulation_functions.R has already been sourced, so that the
# following helpers already exist in the environment:
#   derive_noreturns_hazards()
#   build_states_targets()
#   solve_uh_given_hu()
#   solve_hu_given_uh()
#   loss_one_age()
#   clamp()
#   clamp_t_feasible()
#   project_point_to_line()
#
# Idea:
#   Instead of snapping p0 = (hu0, uh0) to the detected line segment between
#   p1 and p2, use that line only to construct a local search bracket.
#   Then, within that bracket, solve the exact one-age curve numerically:
#     hu -> uh_hat(hu) via solve_uh_given_hu()
#   and minimize squared distance from p0 to the curve point
#     (hu, uh_hat(hu)).
#
# This is intentionally minimal and reuses the existing machinery.

.curve_snap_one_age <- function(i,
                                p0,
                                p1,
                                p2,
                                states,
                                hd,
                                ud,
                                hu_lo,
                                hu_hi,
                                uh_lo,
                                uh_hi,
                                cap,
                                age_int,
                                tiny,
                                curve_tol,
                                curve_pad_mult,
                                curve_pad_abs,
                                curve_loss_w,
                                fallback_to_line = TRUE) {
  
  v <- p2 - p1
  
  # Degenerate / vertical-ish line fallback
  if (!all(is.finite(v)) || sum(v * v) < 1e-30 || abs(v[1]) < 1e-14) {
    if (isTRUE(fallback_to_line)) {
      return(list(
        hu = p1[1],
        uh = p1[2],
        loss = loss_one_age(i, states, hd, ud, hu = p1[1], uh = p1[2], age_int = age_int),
        method = "fallback_p1",
        ok = TRUE
      ))
    }
    return(list(hu = NA_real_, uh = NA_real_, loss = Inf, method = "failed", ok = FALSE))
  }
  
  feas <- clamp_t_feasible(
    a = p1, v = v,
    hu_lo = hu_lo, hu_hi = hu_hi,
    uh_lo = uh_lo, uh_hi = uh_hi,
    cap = cap
  )
  
  if (!isTRUE(feas$ok)) {
    return(list(
      hu = p1[1],
      uh = p1[2],
      loss = loss_one_age(i, states, hd, ud, hu = p1[1], uh = p1[2], age_int = age_int),
      method = "fallback_p1",
      ok = TRUE
    ))
  }
  
  pr <- project_point_to_line(p0, p1, v)
  t_proj <- clamp(pr$t, feas$t_lo, feas$t_hi)
  p_proj <- p1 + t_proj * v
  hu_proj <- p_proj[1]
  
  # Use the feasible line segment only to define a local hu bracket.
  # The actual snap is to the exact curve, not to the line.
  t_end <- c(feas$t_lo, feas$t_hi)
  hu_end <- p1[1] + t_end * v[1]
  hu_seg_lo <- max(hu_lo, min(hu_end), p1[1])
  hu_seg_hi <- min(hu_hi, max(hu_end), cap - uh_lo)
  
  if (!is.finite(hu_seg_lo) || !is.finite(hu_seg_hi) || hu_seg_lo >= hu_seg_hi) {
    hu_seg_lo <- max(hu_lo, p1[1])
    hu_seg_hi <- min(hu_hi, cap - uh_lo)
  }
  
  pad <- max(curve_pad_abs, curve_pad_mult * abs(hu_proj))
  lo <- max(hu_seg_lo, hu_proj - pad)
  hi <- min(hu_seg_hi, hu_proj + pad)
  
  # Ensure bracket is valid; otherwise fall back to the full feasible hu segment.
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    lo <- hu_seg_lo
    hi <- hu_seg_hi
  }
  
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    return(list(
      hu = p_proj[1],
      uh = p_proj[2],
      loss = loss_one_age(i, states, hd, ud, hu = p_proj[1], uh = p_proj[2], age_int = age_int),
      method = "fallback_line",
      ok = TRUE
    ))
  }
  
  obj_hu <- function(hu) {
    uh_hi_eff <- min(uh_hi, cap - hu)
    if (!is.finite(uh_hi_eff) || uh_hi_eff < uh_lo) return(.Machine$double.xmax / 100)
    
    sol <- solve_uh_given_hu(
      i, states, hd, ud,
      hu = hu,
      uh_lo = uh_lo, uh_hi = uh_hi_eff,
      age_int = age_int, tiny = tiny
    )
    
    if (!is.finite(sol$uh) || !is.finite(sol$loss)) return(.Machine$double.xmax / 100)
    
    d2 <- (hu - p0[1])^2 + (sol$uh - p0[2])^2
    d2 + curve_loss_w * sol$loss
  }
  
  opt <- optimize(obj_hu, interval = c(lo, hi), tol = curve_tol)
  hu_star <- unname(opt$minimum)
  uh_hi_eff <- min(uh_hi, cap - hu_star)
  
  sol_star <- solve_uh_given_hu(
    i, states, hd, ud,
    hu = hu_star,
    uh_lo = uh_lo, uh_hi = uh_hi_eff,
    age_int = age_int, tiny = tiny
  )
  
  if (!is.finite(sol_star$uh) || !is.finite(sol_star$loss)) {
    return(list(
      hu = p_proj[1],
      uh = p_proj[2],
      loss = loss_one_age(i, states, hd, ud, hu = p_proj[1], uh = p_proj[2], age_int = age_int),
      method = "fallback_line",
      ok = TRUE
    ))
  }
  
  list(
    hu = hu_star,
    uh = sol_star$uh,
    loss = sol_star$loss,
    method = "curve",
    ok = TRUE,
    hu_proj = hu_proj,
    bracket_lo = lo,
    bracket_hi = hi,
    objective = opt$objective
  )
}


# derive_returns_hazards_from_Rx_curve <- function(
    #     age,
#     lx,
#     prev,
#     Rx,
#     age_int = 1,
#     bounds = c(1e-12, 2),
#     hu_bounds = NULL,
#     uh_bounds = NULL,
#     maxit = 500,
#     reltol = 1e-14,
#     hu_0,
#     uh_0,
#     smooth_w = 0.5,
#     curvature_w = 0.5,
#     bound_w  = 1e-3,
#     turnover_K = 1.2,
#     turnover_w = 1e4,
#     fit_w = 1e6,
#     solver = c("project"),
#     snap_tol = 1e-20,
#     verbose = FALSE,
#     eps_log = 1e-2,
#     tiny = 1e-14,
#     line_tol = 1e-12,
#     final_refine = TRUE,
#     refine_tol = 1e-16,
#     refine_dist_w = 1e-6,
#     extrap_last_age = TRUE,
#     # curve snap controls
#     curve_tol = 1e-10,
#     curve_pad_mult = 0.25,
#     curve_pad_abs = 1e-6,
#     curve_loss_w = 1,
#     fallback_to_line = TRUE
# ){
#   
#   solver <- match.arg(solver)
#   
#   n <- length(age)
#   stopifnot(length(lx) == n, length(prev) == n, length(Rx) == n,
#             length(hu_0) == n, length(uh_0) == n)
#   
#   hz_nr <- derive_noreturns_hazards(
#     age = age,
#     lx = lx,
#     prev = prev,
#     Rx = Rx,
#     age_int = age_int,
#     verbose = FALSE
#   )
#   
#   hd <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
#   ud <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
#   hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
#   
#   if (is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   if (is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   
#   hu_lo <- pmax(hu_bounds[, 1], hu_nr)
#   hu_hi <- hu_bounds[, 2]
#   uh_lo <- pmax(0, uh_bounds[, 1])
#   uh_hi <- uh_bounds[, 2]
#   
#   states <- build_states_targets(age, lx, prev)
#   
#   hu_out <- numeric(n)
#   uh_out <- numeric(n)
#   p2_loss <- rep(NA_real_, n)
#   degenerate <- rep(FALSE, n)
#   final_loss <- rep(NA_real_, n)
#   method <- rep(NA_character_, n)
#   hu_proj_out <- rep(NA_real_, n)
#   bracket_lo <- rep(NA_real_, n)
#   bracket_hi <- rep(NA_real_, n)
#   
#   for (i in seq_len(n)) {
#     
#     cap <- turnover_K * (hu_0[i] + uh_0[i])
#     
#     p1 <- c(clamp(hu_nr[i], hu_lo[i], hu_hi[i]), 0)
#     p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
#     
#     p0 <- c(clamp(hu_0[i], hu_lo[i], hu_hi[i]), clamp(uh_0[i], uh_lo[i], uh_hi[i]))
#     
#     hu2_try <- clamp(p1[1] * exp(eps_log), hu_lo[i], hu_hi[i])
#     uh_hi_eff <- min(uh_hi[i], cap - hu2_try)
#     if (!is.finite(uh_hi_eff)) uh_hi_eff <- uh_hi[i]
#     uh_hi_eff <- max(uh_hi_eff, uh_lo[i])
#     
#     sol2 <- solve_uh_given_hu(
#       i, states, hd, ud,
#       hu = hu2_try,
#       uh_lo = uh_lo[i], uh_hi = uh_hi_eff,
#       age_int = age_int, tiny = tiny
#     )
#     
#     have_p2 <- is.finite(sol2$loss) && is.finite(sol2$uh) &&
#       (abs(hu2_try - p1[1]) > 0) && (sol2$loss <= line_tol)
#     
#     if (!have_p2) {
#       uh_seed <- clamp(max(uh_lo[i], tiny) * exp(eps_log), uh_lo[i], uh_hi[i])
#       hu_hi_eff <- min(hu_hi[i], cap - uh_seed)
#       if (!is.finite(hu_hi_eff)) hu_hi_eff <- hu_hi[i]
#       hu_hi_eff <- max(hu_hi_eff, hu_lo[i])
#       
#       sol2b <- solve_hu_given_uh(
#         i, states, hd, ud,
#         uh = uh_seed,
#         hu_lo = hu_lo[i], hu_hi = hu_hi_eff,
#         age_int = age_int
#       )
#       
#       have_p2 <- is.finite(sol2b$loss) && is.finite(sol2b$hu) &&
#         (abs(uh_seed - p1[2]) > 0) && (sol2b$loss <= line_tol)
#       
#       if (have_p2) {
#         p2 <- c(sol2b$hu, uh_seed)
#         p2_loss[i] <- sol2b$loss
#       } else {
#         if (is.finite(sol2$loss) && is.finite(sol2$uh) && abs(hu2_try - p1[1]) > 0) {
#           p2 <- c(hu2_try, sol2$uh)
#           p2_loss[i] <- sol2$loss
#           have_p2 <- TRUE
#         } else if (is.finite(sol2b$loss) && is.finite(sol2b$hu) && abs(uh_seed - p1[2]) > 0) {
#           p2 <- c(sol2b$hu, uh_seed)
#           p2_loss[i] <- sol2b$loss
#           have_p2 <- TRUE
#         } else {
#           p2 <- p1
#           have_p2 <- FALSE
#         }
#       }
#     } else {
#       p2 <- c(hu2_try, sol2$uh)
#       p2_loss[i] <- sol2$loss
#     }
#     
#     out_i <- .curve_snap_one_age(
#       i = i,
#       p0 = p0,
#       p1 = p1,
#       p2 = p2,
#       states = states,
#       hd = hd,
#       ud = ud,
#       hu_lo = hu_lo[i],
#       hu_hi = hu_hi[i],
#       uh_lo = uh_lo[i],
#       uh_hi = uh_hi[i],
#       cap = cap,
#       age_int = age_int,
#       tiny = tiny,
#       curve_tol = curve_tol,
#       curve_pad_mult = curve_pad_mult,
#       curve_pad_abs = curve_pad_abs,
#       curve_loss_w = curve_loss_w,
#       fallback_to_line = fallback_to_line
#     )
#     
#     hu_out[i] <- out_i$hu
#     uh_out[i] <- out_i$uh
#     final_loss[i] <- out_i$loss
#     method[i] <- out_i$method
#     hu_proj_out[i] <- if (!is.null(out_i$hu_proj)) out_i$hu_proj else NA_real_
#     bracket_lo[i] <- if (!is.null(out_i$bracket_lo)) out_i$bracket_lo else NA_real_
#     bracket_hi[i] <- if (!is.null(out_i$bracket_hi)) out_i$bracket_hi else NA_real_
#     
#     degenerate[i] <- isTRUE(all.equal(p1, p2, tolerance = 1e-15))
#     
#     if (verbose && (i %% 10 == 0 || i == 1 || i == n)) {
#       message(sprintf(
#         "age %s: curve_loss=%.3e  p2_loss=%s  method=%s",
#         age[i], final_loss[i], format(p2_loss[i], scientific = TRUE, digits = 2), method[i]
#       ))
#     }
#   }
#   
#   if (isTRUE(extrap_last_age) && n >= 3) {
#     min_haz <- bounds[1]
#     extrap_last <- function(x) {
#       if (is.finite(x[n - 1]) && is.finite(x[n - 2]) && x[n - 1] > 0 && x[n - 2] > 0) {
#         exp(2 * log(x[n - 1]) - log(x[n - 2]))
#       } else if (is.finite(x[n - 1])) {
#         x[n - 1]
#       } else {
#         min_haz
#       }
#     }
#     
#     hu_out[n] <- clamp(pmax(extrap_last(hu_out), min_haz), hu_lo[n], hu_hi[n])
#     uh_out[n] <- clamp(pmax(extrap_last(uh_out), min_haz), uh_lo[n], uh_hi[n])
#     
#     cap_n <- turnover_K * (hu_0[n] + uh_0[n])
#     if (is.finite(cap_n) && (hu_out[n] + uh_out[n] > cap_n)) {
#       uh_out[n] <- max(uh_lo[n], cap_n - hu_out[n])
#     }
#   }
#   
#   out <- tibble::tibble(
#     age = rep(age, times = 4),
#     trans = rep(c("hd", "hu", "ud", "uh"), each = n),
#     hazard = c(hd, hu_out, ud, uh_out)
#   )
#   
#   attr(out, "diagnostics") <- list(
#     p2_loss = p2_loss,
#     degenerate = degenerate,
#     final_loss = final_loss,
#     method = method,
#     hu_proj = hu_proj_out,
#     bracket_lo = bracket_lo,
#     bracket_hi = bracket_hi
#   )
#   
#   out
# }

# derive_returns_hazards_from_Rx_sr_local_tight.R
#
# Drop-in redux of derive_returns_hazards_from_Rx() for the decomp_sullivan
# simulation workflow.
#
# Assumes simulation_functions.R has already been sourced, so the following
# existing helpers are available and are reused here rather than redefined:
#   - derive_noreturns_hazards()
#   - build_states_targets()
#   - solve_uh_given_hu()
#   - solve_hu_given_uh()
#   - loss_one_age()
#   - project_point_to_line()
#   - clamp_t_feasible()
#   - clamp()
#
# Main change versus the earlier sr_local_tight version:
#   the local s/r objective is written explicitly in terms of the exact
#   next-age lH/lU target rather than calling loss_one_age().
#   This is mathematically equivalent to loss_one_age() in
#   simulation_functions.R and is included only to make the target explicit.

.sr_to_haz <- function(par) {
  s <- par[1]
  r <- par[2]
  T <- exp(s)
  q <- exp(r)
  uh <- T / (1 + q)
  hu <- T - uh
  c(hu = hu, uh = uh)
}

.sr_objective_local <- function(par,
                                i, states, hd, ud,
                                hu_lo, hu_hi,
                                uh_lo, uh_hi,
                                cap,
                                age_int,
                                par0 = NULL,
                                anchor_w = 0,
                                obj_scale = 1e12) {
  hz <- .sr_to_haz(par)
  hu <- hz["hu"]
  uh <- hz["uh"]
  
  if (!is.finite(hu) || !is.finite(uh) || hu <= 0 || uh < 0) {
    return(1e100)
  }
  
  pen <- 0
  if (hu < hu_lo) pen <- pen + (hu_lo - hu)^2
  if (hu > hu_hi) pen <- pen + (hu - hu_hi)^2
  if (uh < uh_lo) pen <- pen + (uh_lo - uh)^2
  if (uh > uh_hi) pen <- pen + (uh - uh_hi)^2
  if (is.finite(cap) && (hu + uh > cap)) pen <- pen + (hu + uh - cap)^2
  
  # Exact one-age target: match next-age lH and lU directly.
  # Note: this is mathematically equivalent to loss_one_age() in
  # simulation_functions.R; it is written out here only for clarity.
  cur <- c(states$lH[i], states$lU[i], states$lD[i])
  targ <- c(states$nextHU[i, 1], states$nextHU[i, 2])
  Q <- make_Q_returns(hu = hu, hd = hd[i], uh = uh, ud = ud[i])
  P <- expm::expm(Q * age_int)
  nxt <- as.numeric(cur %*% P)
  raw_obj <- sum((nxt[1:2] - targ)^2)
  
  if (pen > 0) {
    raw_obj <- raw_obj + 1e12 * pen
  }
  
  if (!is.null(par0) && is.finite(anchor_w) && anchor_w > 0) {
    raw_obj <- raw_obj + anchor_w * sum((par - par0)^2)
  }
  
  obj_scale * raw_obj
}

# derive_returns_hazards_from_Rx <- function(
    #     age,
#     lx,
#     prev,
#     Rx,
#     age_int = 1,
#     bounds = c(1e-12, 2),
#     hu_bounds = NULL,
#     uh_bounds = NULL,
#     maxit = 500,
#     reltol = 1e-14,
#     hu_0,
#     uh_0,
#     smooth_w = 0.5,
#     curvature_w = 0.5,
#     bound_w  = 1e-3,
#     turnover_K = 1.2,
#     turnover_w = 1e4,
#     fit_w = 1e6,
#     solver = c("project"),
#     snap_tol = 1e-20,
#     verbose = FALSE,
#     eps_log = 1e-2,
#     tiny = 1e-14,
#     line_tol = 1e-12,
#     final_refine = TRUE,
#     refine_tol = 1e-16,
#     refine_dist_w = 1e-6,
#     extrap_last_age = TRUE,
#     # local s/r refinement controls
#     sr_refine = TRUE,
#     sr_skip_last_age = TRUE,
#     sr_improve_tol = 0,
#     sr_s_window = log(1.25),
#     sr_r_window = log(1.5),
#     sr_anchor_w = 0,
#     sr_maxit = 200,
#     sr_factr = 1,
#     sr_pgtol = 0,
#     sr_obj_scale = 1e12
# ){
#   solver <- match.arg(solver)
#   
#   age  <- unname(as.numeric(age))
#   lx   <- unname(as.numeric(lx))
#   prev <- unname(as.numeric(prev))
#   Rx   <- unname(as.numeric(Rx))
#   hu_0 <- unname(as.numeric(hu_0))
#   uh_0 <- unname(as.numeric(uh_0))
#   
#   n <- length(age)
#   stopifnot(length(lx) == n, length(prev) == n, length(Rx) == n,
#             length(hu_0) == n, length(uh_0) == n)
#   
#   hz_nr <- derive_noreturns_hazards(
#     age = age,
#     lx = lx,
#     prev = prev,
#     Rx = Rx,
#     age_int = age_int,
#     verbose = FALSE
#   )
#   
#   hd <- hz_nr[hz_nr$trans == "hd", "hazard", drop = TRUE]
#   ud <- hz_nr[hz_nr$trans == "ud", "hazard", drop = TRUE]
#   hu_nr <- hz_nr[hz_nr$trans == "hu", "hazard", drop = TRUE]
#   
#   if (is.null(hu_bounds)) hu_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   if (is.null(uh_bounds)) uh_bounds <- cbind(rep(bounds[1], n), rep(bounds[2], n))
#   
#   hu_lo <- pmax(hu_bounds[,1], hu_nr)
#   hu_hi <- hu_bounds[,2]
#   uh_lo <- pmax(0, uh_bounds[,1])
#   uh_hi <- uh_bounds[,2]
#   
#   states <- build_states_targets(age, lx, prev)
#   
#   hu_out <- numeric(n)
#   uh_out <- numeric(n)
#   p2_loss <- rep(NA_real_, n)
#   degenerate <- rep(FALSE, n)
#   pre_sr_loss <- rep(NA_real_, n)
#   post_sr_loss <- rep(NA_real_, n)
#   accepted_sr <- rep(FALSE, n)
#   sr_counts <- rep(NA_integer_, n)
#   sr_convergence <- rep(NA_integer_, n)
#   sr_message <- rep(NA_character_, n)
#   sr_hu <- rep(NA_real_, n)
#   sr_uh <- rep(NA_real_, n)
#   
#   for (i in seq_len(n)) {
#     cap <- turnover_K * (hu_0[i] + uh_0[i])
#     
#     p1 <- c(clamp(hu_nr[i], hu_lo[i], hu_hi[i]), 0)
#     p1[2] <- clamp(p1[2], uh_lo[i], uh_hi[i])
#     
#     p0 <- c(clamp(hu_0[i], hu_lo[i], hu_hi[i]), clamp(uh_0[i], uh_lo[i], uh_hi[i]))
#     
#     hu2_try <- clamp(p1[1] * exp(eps_log), hu_lo[i], hu_hi[i])
#     uh_hi_eff <- min(uh_hi[i], cap - hu2_try)
#     if (!is.finite(uh_hi_eff)) uh_hi_eff <- uh_hi[i]
#     uh_hi_eff <- max(uh_hi_eff, uh_lo[i])
#     
#     sol2 <- solve_uh_given_hu(
#       i, states, hd, ud,
#       hu = hu2_try,
#       uh_lo = uh_lo[i], uh_hi = uh_hi_eff,
#       age_int = age_int, tiny = tiny
#     )
#     
#     have_p2 <- is.finite(sol2$loss) && is.finite(sol2$uh) &&
#       (abs(hu2_try - p1[1]) > 0) && (sol2$loss <= line_tol)
#     
#     if (!have_p2) {
#       uh_seed <- clamp(max(uh_lo[i], tiny) * exp(eps_log), uh_lo[i], uh_hi[i])
#       hu_hi_eff <- min(hu_hi[i], cap - uh_seed)
#       if (!is.finite(hu_hi_eff)) hu_hi_eff <- hu_hi[i]
#       hu_hi_eff <- max(hu_hi_eff, hu_lo[i])
#       
#       sol2b <- solve_hu_given_uh(
#         i, states, hd, ud,
#         uh = uh_seed,
#         hu_lo = hu_lo[i], hu_hi = hu_hi_eff,
#         age_int = age_int
#       )
#       
#       have_p2 <- is.finite(sol2b$loss) && is.finite(sol2b$hu) &&
#         (abs(uh_seed - p1[2]) > 0) && (sol2b$loss <= line_tol)
#       
#       if (have_p2) {
#         p2 <- c(sol2b$hu, uh_seed)
#         p2_loss[i] <- sol2b$loss
#       } else {
#         if (is.finite(sol2$loss) && is.finite(sol2$uh) && abs(hu2_try - p1[1]) > 0) {
#           p2 <- c(hu2_try, sol2$uh)
#           p2_loss[i] <- sol2$loss
#           have_p2 <- TRUE
#         } else if (is.finite(sol2b$loss) && is.finite(sol2b$hu) && abs(uh_seed - p1[2]) > 0) {
#           p2 <- c(sol2b$hu, uh_seed)
#           p2_loss[i] <- sol2b$loss
#           have_p2 <- TRUE
#         } else {
#           p2 <- p1
#           have_p2 <- FALSE
#         }
#       }
#     } else {
#       p2 <- c(hu2_try, sol2$uh)
#       p2_loss[i] <- sol2$loss
#     }
#     
#     v <- p2 - p1
#     if (sum(v*v) < 1e-30) {
#       degenerate[i] <- TRUE
#       hu_hat <- p1[1]
#       uh_hat <- p1[2]
#     } else {
#       feas <- clamp_t_feasible(
#         a = p1, v = v,
#         hu_lo = hu_lo[i], hu_hi = hu_hi[i],
#         uh_lo = uh_lo[i], uh_hi = uh_hi[i],
#         cap = cap
#       )
#       
#       if (!isTRUE(feas$ok)) {
#         hu_hat <- p1[1]
#         uh_hat <- p1[2]
#       } else {
#         pr <- project_point_to_line(p0, p1, v)
#         t_proj <- clamp(pr$t, feas$t_lo, feas$t_hi)
#         
#         if (isTRUE(final_refine)) {
#           f_t <- function(t) {
#             pt <- p1 + t * v
#             ls <- loss_one_age(i, states, hd, ud, hu = pt[1], uh = pt[2], age_int = age_int)
#             ls + refine_dist_w * (t - t_proj)^2
#           }
#           opt_t <- optimize(f_t, interval = c(feas$t_lo, feas$t_hi), tol = refine_tol)
#           t_star <- opt_t$minimum
#         } else {
#           t_star <- t_proj
#         }
#         
#         p_star <- p1 + t_star * v
#         hu_hat <- p_star[1]
#         uh_hat <- p_star[2]
#       }
#     }
#     
#     pre_sr_loss[i] <- loss_one_age(i, states, hd, ud, hu = hu_hat, uh = uh_hat, age_int = age_int)
#     
#     do_sr <- isTRUE(sr_refine) && !(isTRUE(sr_skip_last_age) && i == n) &&
#       is.finite(hu_hat) && is.finite(uh_hat) && hu_hat > 0 && uh_hat > 0
#     
#     if (do_sr) {
#       s0 <- log(hu_hat + uh_hat)
#       r0 <- log(hu_hat / uh_hat)
#       par0 <- c(s0, r0)
#       
#       lower <- c(s0 - sr_s_window, r0 - sr_r_window)
#       upper <- c(s0 + sr_s_window, r0 + sr_r_window)
#       
#       obj_sr <- function(par) {
#         .sr_objective_local(
#           par = par,
#           i = i, states = states, hd = hd, ud = ud,
#           hu_lo = hu_lo[i], hu_hi = hu_hi[i],
#           uh_lo = uh_lo[i], uh_hi = uh_hi[i],
#           cap = cap,
#           age_int = age_int,
#           par0 = par0,
#           anchor_w = sr_anchor_w,
#           obj_scale = sr_obj_scale
#         )
#       }
#       
#       opt_sr <- optim(
#         par = par0,
#         fn = obj_sr,
#         method = "L-BFGS-B",
#         lower = lower,
#         upper = upper,
#         control = list(
#           maxit = sr_maxit,
#           factr = sr_factr,
#           pgtol = sr_pgtol,
#           trace = 0,
#           REPORT = 10
#         )
#       )
#       
#       sr_counts[i] <- opt_sr$counts[[1]]
#       sr_convergence[i] <- opt_sr$convergence
#       if (!is.null(opt_sr$message)) sr_message[i] <- as.character(opt_sr$message)[1]
#       
#       hz_sr <- .sr_to_haz(opt_sr$par)
#       hu_try <- as.numeric(hz_sr["hu"])
#       uh_try <- as.numeric(hz_sr["uh"])
#       loss_try <- loss_one_age(i, states, hd, ud, hu = hu_try, uh = uh_try, age_int = age_int)
#       
#       sr_hu[i] <- hu_try
#       sr_uh[i] <- uh_try
#       post_sr_loss[i] <- loss_try
#       
#       if (is.finite(loss_try) && (loss_try + sr_improve_tol < pre_sr_loss[i])) {
#         hu_hat <- hu_try
#         uh_hat <- uh_try
#         accepted_sr[i] <- TRUE
#       }
#     } else {
#       post_sr_loss[i] <- pre_sr_loss[i]
#     }
#     
#     hu_out[i] <- hu_hat
#     uh_out[i] <- uh_hat
#     
#     if (verbose && (i %% 10 == 0 || i == 1 || i == n)) {
#       ls <- loss_one_age(i, states, hd, ud, hu = hu_out[i], uh = uh_out[i], age_int = age_int)
#       message(sprintf(
#         "age %s: loss=%.3e  pre_sr=%.3e  post_sr=%s  p2_loss=%s  acc=%s  deg=%s",
#         age[i], ls, pre_sr_loss[i],
#         format(post_sr_loss[i], scientific = TRUE, digits = 2),
#         format(p2_loss[i], scientific = TRUE, digits = 2),
#         accepted_sr[i], degenerate[i]
#       ))
#     }
#   }
#   
#   if (isTRUE(extrap_last_age) && n >= 3) {
#     min_haz <- bounds[1]
#     extrap_last <- function(x) {
#       if (is.finite(x[n-1]) && is.finite(x[n-2]) && x[n-1] > 0 && x[n-2] > 0) {
#         exp(2 * log(x[n-1]) - log(x[n-2]))
#       } else if (is.finite(x[n-1])) {
#         x[n-1]
#       } else {
#         min_haz
#       }
#     }
#     
#     hu_out[n] <- clamp(pmax(extrap_last(hu_out), min_haz), hu_lo[n], hu_hi[n])
#     uh_out[n] <- clamp(pmax(extrap_last(uh_out), min_haz), uh_lo[n], uh_hi[n])
#     
#     cap_n <- turnover_K * (hu_0[n] + uh_0[n])
#     if (is.finite(cap_n) && (hu_out[n] + uh_out[n] > cap_n)) {
#       uh_out[n] <- max(uh_lo[n], cap_n - hu_out[n])
#     }
#   }
#   
#   out <- tibble::tibble(
#     age = rep(age, times = 4),
#     trans = rep(c("hd", "hu", "ud", "uh"), each = n),
#     hazard = c(hd, hu_out, ud, uh_out)
#   )
#   
#   attr(out, "diagnostics") <- list(
#     p2_loss = p2_loss,
#     degenerate = degenerate,
#     pre_sr_loss = pre_sr_loss,
#     post_sr_loss = post_sr_loss,
#     accepted_sr = accepted_sr,
#     sr_counts = sr_counts,
#     sr_convergence = sr_convergence,
#     sr_message = sr_message,
#     sr_hu = sr_hu,
#     sr_uh = sr_uh
#   )
#   
#   out
# }