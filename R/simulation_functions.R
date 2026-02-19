# ------------------------------------------------------------------------------
# simulation_functions2.R (redux)
#
# Create multiple multistate "worlds" (hazard schedules) that reproduce nearly
# identical Sullivan outputs (lx and prevalence), via:
#   (1) Basic functionality (hazards -> probs -> multistate lifetable)
#   (2) SVD-based near-null perturbations (local non-identifiability directions)
#   (3) Final polishing via ridge-regularized optimization with sign constraints
#       and optional, tame piecewise-loglinear "kink" at a fixed knot age.
#
# Parameterizations
#   Linear (legacy): 8 params
#     c(i_hu,i_hd,i_uh,i_ud, s_hu,s_hd,s_uh,s_ud)
#
#   Piecewise log-linear with 1 knot (recommended for extra flexibility): 12 params
#     c(i_hu,i_hd,i_uh,i_ud,
#       s1_hu,s1_hd,s1_uh,s1_ud,
#       s2_hu,s2_hd,s2_uh,s2_ud)
#     where s1 = slope before knot, s2 = slope after knot.
#     Setting s2_* == s1_* gives exactly linear hazards (no kink).
#
# Knot age is fixed (default 75), same for all transitions.
# ------------------------------------------------------------------------------

# ---------
# (1) Basic functionality: parameter helpers + hazards -> probs -> lifetable
# ---------

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

# Accept either 8 (linear) or 12 (piecewise) parameter vectors.
# Returns a named vector in canonical order.
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
  
  # If mismatch (e.g., asked piecewise but got 8), just return canonical for its length.
  if (length(par_vec) == 8) return(par_vec[ms_par_names(FALSE)])
  par_vec[ms_par_names(TRUE)]
}

# Convert 8->12 by duplicating slopes: s1 = s, s2 = s
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

# Convert 12->8 by keeping s_hu = s1_hu (and dropping s2)
# (used rarely; mostly you stay in 12 once you add kinks)
as_linear_pars <- function(pars) {
  v <- validate_par_vec(as_par_vec(pars), piecewise = NULL)
  if (length(v) == 8) return(v)
  
  out <- c(
    v[c("i_hu","i_hd","i_uh","i_ud")],
    setNames(as.numeric(v[c("s1_hu","s1_hd","s1_uh","s1_ud")]), paste0("s_", c("hu","hd","uh","ud")))
  )
  out[ms_par_names(FALSE)]
}

pars_vec_to_df <- function(par_vec) {
  v <- validate_par_vec(par_vec, piecewise = NULL)
  if (length(v) == 8) {
    tibble::tibble(
      trans     = c("hu","hd","uh","ud"),
      intercept = as.numeric(v[c("i_hu","i_hd","i_uh","i_ud")]),
      slope     = as.numeric(v[c("s_hu","s_hd","s_uh","s_ud")])
    )
  } else {
    tibble::tibble(
      trans     = c("hu","hd","uh","ud"),
      intercept = as.numeric(v[c("i_hu","i_hd","i_uh","i_ud")]),
      slope1    = as.numeric(v[c("s1_hu","s1_hd","s1_uh","s1_ud")]),
      slope2    = as.numeric(v[c("s2_hu","s2_hd","s2_uh","s2_ud")])
    )
  }
}

# Allow df inputs too (legacy). For df, assume linear unless slope1/slope2 present.
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

as_par_df <- function(pars) pars_vec_to_df(as_par_vec(pars))

# Expand hazards by age.
# Linear: log h = i + s * x
# Piecewise: log h = i + s1 * x + (s2 - s1) * max(0, age_mid - knot_age)
# where x = age_mid - min(age), age_mid = age + age_int/2
expand_hazards <- function(pars,
                           age = 50:100,
                           age_int = 1,
                           knot_age = 75) {
  
  v <- as_par_vec(pars)
  age_mid <- age + age_int / 2
  x <- age_mid - min(age)
  
  if (length(v) == 8) {
    pars_df <- pars_vec_to_df(v)
    tidyr::expand_grid(age = age, trans = pars_df$trans) |>
      dplyr::left_join(pars_df, by = "trans") |>
      dplyr::mutate(
        log_hazard = intercept + slope * (age + age_int/2 - min(age)),
        hazard = exp(log_hazard)
      ) |>
      dplyr::select(age, trans, hazard)
  } else {
    # piecewise
    i <- v[c("i_hu","i_hd","i_uh","i_ud")]
    s1 <- v[c("s1_hu","s1_hd","s1_uh","s1_ud")]
    s2 <- v[c("s2_hu","s2_hd","s2_uh","s2_ud")]
    
    trans <- c("hu","hd","uh","ud")
    i <- setNames(as.numeric(i), trans)
    s1 <- setNames(as.numeric(s1), trans)
    s2 <- setNames(as.numeric(s2), trans)
    
    kink <- pmax(0, age_mid - knot_age)
    
    tidyr::expand_grid(age = age, trans = trans) |>
      dplyr::mutate(
        age_mid = age + age_int/2,
        x = age_mid - min(age),
        kink = pmax(0, age_mid - knot_age),
        intercept = i[trans],
        slope1 = s1[trans],
        slope2 = s2[trans],
        log_hazard = intercept + slope1 * x + (slope2 - slope1) * kink,
        hazard = exp(log_hazard)
      ) |>
      dplyr::select(age, trans, hazard)
  }
}

# Build CTMC generator Q from hazards (hu, hd, uh, ud)
make_Q <- function(x) {
  with(as.list(x),
       matrix(c(
         -(hu + hd),  hu,  hd,
         uh, -(uh + ud), ud,
         0,   0,   0
       ), nrow = 3, byrow = TRUE)
  )
}

# Convert hazards to transition probability matrices via P = expm(Q * age_int)
haz_to_probs <- function(hazard, age, age_int = 1) {
  if (length(age_int) != 1 || !is.numeric(age_int) || age_int <= 0) {
    stop("age_int must be a single positive numeric value.", call. = FALSE)
  }
  
  hazard |>
    dplyr::group_nest(age) |>
    dplyr::mutate(
      Q = purrr::map(.data$data, ~ make_Q(tibble::deframe(.x))),
      P = purrr::map(.data$Q, ~ (.x * age_int) |>
                       expm::expm() |>
                       as.data.frame() |>
                       rlang::set_names(c("H","U","D")) |>
                       cbind("from" = c("H","U","D")))
    ) |>
    dplyr::select(age, P) |>
    tidyr::unnest(P)
}

# Multistate lifetable recursion (init defaults to c(H=1,U=0))
calculate_lt <- function(P, init = c(H = 1, U = 0)) {
  P |>
    tidyr::pivot_longer(-c(age, from), names_to = "to", values_to = "p") |>
    dplyr::filter(.data$from != "D") |>
    tidyr::unite("state", c("from","to"), sep = "") |>
    tidyr::pivot_wider(names_from = state, values_from = p) |>
    dplyr::group_modify(~{
      .x$lh <- numeric(nrow(.x))
      .x$lu <- numeric(nrow(.x))
      
      init_use <- init
      .x$lh[1] <- init_use["H"]
      .x$lu[1] <- init_use["U"]
      
      for (ii in seq_len(nrow(.x) - 1)) {
        .x$lh[ii + 1] <- .x$lh[ii] * .x$HH[ii] + .x$lu[ii] * .x$UH[ii]
        .x$lu[ii + 1] <- .x$lu[ii] * .x$UU[ii] + .x$lh[ii] * .x$HU[ii]
      }
      .x
    }) |>
    dplyr::mutate(
      prevalence = .data$lu / (.data$lh + .data$lu),
      lx         = .data$lh + .data$lu,
      Lx         = (.data$lx + dplyr::lead(.data$lx, default = 0)) / 2
    )
}

run_mslt <- function(pars,
                     age,
                     init = c(H = 1, U = 0),
                     age_int = 1,
                     knot_age = 75) {
  haz <- expand_hazards(pars, age = age, age_int = age_int, knot_age = knot_age)
  P   <- haz_to_probs(haz, age = age, age_int = age_int)
  calculate_lt(P, init = init)
}

# ---------
# (2) SVD-based perturbations: find near-null directions and apply them
# ---------

mslt_summary_vector <- function(pars,
                                age,
                                init = c(H = 1, U = 0),
                                age_int = 1,
                                knot_age = 75,
                                outputs = c("lx", "prevalence"),
                                weights = NULL,
                                standardize = TRUE) {
  
  lt <- run_mslt(pars, age = age, init = init, age_int = age_int, knot_age = knot_age)
  
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

# A Jacobian  is a matrix composed of all first-order partial derivatives 
# of a vector-valued function, mapping an n-dimensional input to an 
# n-dimensional output. It represents the best linear approximation of a 
# nonlinear function at a specific point, acting as a multivariable
# generalization of the derivative. The "central" direction is basically
# doing the same operation as inside the horiuchi method, except at just 
# one point, and for many outputs.
mslt_jacobian <- function(base_pars,
                          free_names,
                          age,
                          init = c(H = 1, U = 0),
                          age_int = 1,
                          knot_age = 75,
                          outputs = c("lx", "prevalence"),
                          weights = NULL,
                          standardize = TRUE,
                          eps = 1e-4,
                          method = c("central", "forward"),
                          step_scale = c("relative", "absolute")) {
  
  method     <- match.arg(method)
  step_scale <- match.arg(step_scale)
  
  theta0 <- as_par_vec(base_pars)
  
  y0 <- mslt_summary_vector(theta0, age, init, age_int, knot_age, outputs, weights, standardize)
  m  <- length(y0)
  p  <- length(free_names)
  
  J <- matrix(NA_real_, nrow = m, ncol = p, dimnames = list(NULL, free_names))
  
  for (j in seq_len(p)) {
    nm <- free_names[j]
    
    step <- if (step_scale == "relative") eps * (abs(theta0[[nm]]) + 1) else eps
    
    theta_plus <- theta0
    theta_plus[[nm]] <- theta0[[nm]] + step
    y_plus <- mslt_summary_vector(theta_plus, age, init, age_int, knot_age, outputs, weights, standardize)
    
    if (method == "central") {
      theta_minus <- theta0
      theta_minus[[nm]] <- theta0[[nm]] - step
      y_minus <- mslt_summary_vector(theta_minus, age, init, age_int, knot_age, outputs, weights, standardize)
      J[, j] <- (y_plus - y_minus) / (2 * step)
    } else {
      J[, j] <- (y_plus - y0) / step
    }
  }
  
  list(J = J, y0 = y0, theta0 = theta0[free_names], free_names = free_names,
       meta = list(age=age, age_int=age_int, knot_age=knot_age, outputs=outputs,
                   standardize=standardize, eps=eps, method=method, step_scale=step_scale))
}
# ---------------------------------------------------------------------------
# This is step 1 of making worlds. We produce pretty darn good worlds here
# and then in a later step we refine (polish) them to make lx and prev match
# even better.
#
# Locally, ms->sullivan mapping can be approximated by its Jacobian:
#     J = d( summaries ) / d( parameters )
#
# If a parameter perturbation vector v satisfies
#     J %*% v ≈ 0
# then moving parameters slightly along v leaves the summaries almost
# unchanged to first order.  Such vectors give *near-null directions*:
# directions in parameter space along which we can move to create alternative
# ms worlds with nearly identical lx and prevalence.
#
# We compute these directions using an SVD of the Jacobian:
#     J = U D V'
# The right singular vectors (columns of V) associated with the smallest
# singular values correspond to directions where the summaries change least.
#
# These vectors provide systematic perturbations of the baseline hazards,
# generating alternative multistate structures that reproduce the same
# aggregate lifetable and prevalence patterns. Or almost anyway.
# ---------------------------------------------------------------------------
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
                                 knot_age = 75,
                                 outputs = c("lx","prevalence"),
                                 tol = 1e-3,
                                 eps = 1e-4,
                                 method = c("central","forward"),
                                 step_scale = c("relative","absolute"),
                                 standardize = TRUE) {
  
  jac <- mslt_jacobian(
    base_pars = base_pars,
    free_names = free_names,
    age = age,
    init = init,
    age_int = age_int,
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

# ---------
# (3) Final polishing: ridge-regularized optimization with sign constraints
# ---------

# Bounds for slopes. For piecewise, enforce sign on BOTH s1 and s2.
# For linear, enforce sign on s.
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

# Objective: Sullivan fit + ridge to stay near theta_target + (optional) 
# kink penalty. kink penalty acts on (s2 - s1) for each transition, keeping
# curvature minor. We penalize distance from the starting parameters
mslt_objective_ridge <- function(par_free,
                                 free_names,
                                 base_pars,
                                 ground_summary,
                                 age,
                                 init = c(H=1,U=0),
                                 age_int = 1,
                                 knot_age = 75,
                                 fixed = NULL,
                                 w_lx = 1,
                                 w_prev = 1,
                                 lambda = 0.3,
                                 theta_target = NULL,
                                 par_scale = NULL,
                                 lambda_kink = 0,
                                 kink_trans = c("hu","hd","uh","ud")) {
  
  if (length(par_free) != length(free_names)) stop("par_free and free_names must have same length.", call. = FALSE)
  
  cand <- as_par_vec(base_pars)
  
  if (!is.null(fixed)) {
    if (is.null(names(fixed)) || any(names(fixed) == "")) stop("fixed must be named.", call. = FALSE)
    cand[names(fixed)] <- as.numeric(fixed)
  }
  
  cand[free_names] <- as.numeric(par_free)
  
  lt_mod <- run_mslt(cand, age = age, init = init, age_int = age_int, knot_age = knot_age) |>
    dplyr::select(age, lx1 = lx, prevalence1 = prevalence) |>
    dplyr::left_join(ground_summary, by = "age") |>
    dplyr::mutate(
      loss_lx   = (lx1 - lx)^2,
      loss_prev = (prevalence1 - prevalence)^2
    )
  
  fit_loss <- sum(w_lx * lt_mod$loss_lx + w_prev * lt_mod$loss_prev)
  
  # ridge identity term
  reg_loss <- 0
  if (!is.null(lambda) && lambda > 0) {
    if (is.null(theta_target)) theta_target <- as_par_vec(base_pars)
    theta_target <- as_par_vec(theta_target)
    
    if (is.null(par_scale)) {
      # robust-ish scaling: per-parameter absolute scale
      par_scale <- rep(1, length(free_names))
      names(par_scale) <- free_names
      for (nm in free_names) {
        par_scale[nm] <- max(1e-6, abs(theta_target[[nm]]) + 1)
      }
    } else {
      if (is.null(names(par_scale))) names(par_scale) <- free_names
      par_scale <- par_scale[free_names]
      par_scale[!is.finite(par_scale) | par_scale <= 0] <- 1
    }
    
    diff <- (cand[free_names] - theta_target[free_names]) / par_scale
    reg_loss <- lambda * sum(as.numeric(diff)^2)
  }
  
  # kink penalty: encourage s2 == s1
  kink_loss <- 0
  if (!is.null(lambda_kink) && lambda_kink > 0) {
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

# Polish a candidate using L-BFGS-B bounds + ridge
# This can take a long time to run since we're optimizing
# many parameters.
polish_candidate_ridge <- function(theta_start,
                                   ground_summary,
                                   free_polish,
                                   age,
                                   init = c(H=1,U=0),
                                   age_int = 1,
                                   knot_age = 75,
                                   fixed = NULL,
                                   lambda = 0.3,
                                   theta_target = theta_start,
                                   lambda_kink = 0,
                                   slope_floor = 1e-6,
                                   maxit = 800) {
  
  b <- make_bounds(free_polish, slope_floor = slope_floor)
  
  res <- optim(
    par = as.numeric(theta_start[free_polish]),
    fn  = mslt_objective_ridge,
    method = "L-BFGS-B",
    lower = b$lower,
    upper = b$upper,
    control = list(maxit = maxit),
    free_names = free_polish,
    base_pars = theta_start,
    ground_summary = ground_summary,
    age = age,
    init = init,
    age_int = age_int,
    knot_age = knot_age,
    fixed = fixed,
    lambda = lambda,
    theta_target = theta_target,
    lambda_kink = lambda_kink
  )
  
  # reconstruct
  cand <- as_par_vec(theta_start)
  if (!is.null(fixed)) cand[names(fixed)] <- as.numeric(fixed)
  cand[free_polish] <- res$par
  cand
}

# ---------
# Utilities: scoring + diversity (used in simulation.R)
# ---------

score_fit <- function(theta, ground_summary, age, init=c(H=1,U=0), age_int=1, knot_age=75) {
  lt <- run_mslt(theta, age=age, init=init, age_int=age_int, knot_age=knot_age)
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

haz_matrix <- function(theta, age, age_int=1, knot_age=75) {
  expand_hazards(theta, age=age, age_int=age_int, knot_age=knot_age) |>
    dplyr::mutate(loghaz = log(hazard)) |>
    dplyr::arrange(trans, age) |>
    dplyr::pull(loghaz)
}

select_farthest <- function(D, k, start = 1) {
  n <- nrow(D)
  sel <- c(start)
  while (length(sel) < k) {
    mind <- apply(D[, sel, drop=FALSE], 1, min)
    mind[sel] <- -Inf
    sel <- c(sel, which.max(mind))
  }
  sel
}