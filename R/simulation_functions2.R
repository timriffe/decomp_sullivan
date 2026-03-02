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
        kink_term = dplyr::case_when(
          haz_shape == "piecewise" ~ pmax(0, age_mid - pivot_age),
          haz_shape == "softkink"  ~ shape_width * (softplus((age_mid - pivot_age)/shape_width) - softplus(z0)),
          TRUE                     ~ 0
        ),
        intercept = i[trans],
        slope1    = s1[trans],
        slope2    = s2[trans],
        log_hazard = intercept + slope1 * x + (slope2 - slope1) * kink_term,
        hazard     = exp(log_hazard)
      ) |>
      dplyr::select(age, trans, hazard)
  }
}
# Q from hazards (hu, hd, uh, ud)
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
                     haz_shape = c("loglinear","piecewise","softkink"),
                     pivot_age = 75,
                     shape_width = 5,
                     knot_age = NULL) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  haz <- expand_hazards(pars, age = age, age_int = age_int, haz_shape = haz_shape,
                        pivot_age = pivot_age, shape_width = shape_width)
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

# Polish a candidate using L-BFGS-B bounds + ridge
# This can take a long time to run since we're optimizing
# many parameters.
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
                                   maxit = 1500,
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
  
  theta1 <- run_one(theta_start, theta_target, lambda, maxit,
                    lambda_kink1, kink_trans1, factr1)
  
  if (!is.null(lambda2) && is.finite(lambda2) && lambda2 >= 0) {
    theta2 <- run_one(theta1, theta1, lambda2, maxit2,
                      lambda_kink2, kink_trans2, factr2)
    return(theta2)
  }
  
  theta1
}
# ---------
# Utilities: scoring + diversity (used in simulation.R)
# ---------

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

haz_matrix <- function(theta, age, age_int = 1,
                       haz_shape = c("loglinear","piecewise","softkink"),
                       pivot_age = 75,
                       shape_width = 5,
                       knot_age = NULL) {
  if (!is.null(knot_age)) pivot_age <- knot_age
  haz_shape <- match.arg(haz_shape)
  expand_hazards(theta, age = age, age_int = age_int, haz_shape = haz_shape,
                 pivot_age = pivot_age, shape_width = shape_width) |>
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

# redux of previous work to derive no-returns worlds
# using an Rx (or Ra) assumption.
# -----------------------------------------------------------------------------
# --------- (3) Derive no-returns hazards from lx, prevalence, and Rx ----------
# Discrete-time (interval) solver consistent with P = expm(Q * age_int).
#
# Goal: Given target lx(x) and prevalence pi(x), plus an assumed hazard ratio
#       Rx(x) = ud(x) / hd(x), construct a no-returns CTMC (uh = 0) such that
#       the one-step update using expm matches BOTH:
#         lx(x+dt) and lu(x+dt) = pi(x+dt)*lx(x+dt)
#
# Unknowns per age interval: hd(x) and hu(x) (with ud(x)=Rx(x)*hd(x), uh(x)=0).
# We solve using nested 1D root-finds:
#   - inner: choose hd to match total survivorship at x+dt
#   - outer: choose hu so that (with that hd) we match lu at x+dt
#
# This avoids prevalence-derivative approximations and typically yields
# essentially exact fits (up to root tolerance).
# -----------------------------------------------------------------------------
derive_noreturns_hazards <- function(age,
                                     lx,
                                     prev,
                                     Rx,
                                     age_int = 1,
                                     min_haz = 1e-12,
                                     # root finding controls
                                     tol = 1e-10,
                                     maxiter = 100,
                                     # brackets (expanded automatically if needed)
                                     hu_bracket = c(1e-12, 5),
                                     hd_bracket = c(1e-12, 5),
                                     verbose = TRUE) {
  
  age <- as.numeric(age)
  n <- length(age)
  if (n < 2) stop("Need at least 2 ages.", call. = FALSE)
  
  lx <- as.numeric(lx)
  if (length(lx) != n) stop("lx must have same length as age.", call. = FALSE)
  
  prev <- as.numeric(prev)
  if (length(prev) != n) stop("prev must have same length as age.", call. = FALSE)
  
  # Rx can be scalar or vector
  if (length(Rx) == 1) Rx <- rep(as.numeric(Rx), n)
  Rx <- as.numeric(Rx)
  if (length(Rx) != n) stop("Rx must be length 1 or length(age).", call. = FALSE)
  
  # target state occupancies
  lu_tgt <- prev * lx
  lh_tgt <- lx - lu_tgt
  
  # numeric safety
  Rx <- pmax(Rx, min_haz)
  lx <- pmax(lx, min_haz)
  lu_tgt <- pmin(pmax(lu_tgt, 0), lx)
  lh_tgt <- pmax(lh_tgt, 0)
  
  # helper: build P given hu, hd, Rx at this age
  make_P <- function(hu, hd, Rx_i, dt) {
    ud <- Rx_i * hd
    uh <- min_haz # effectively zero returns
    Q <- matrix(c(
      -(hu + hd),  hu,  hd,
      uh, -(uh + ud), ud,
      0, 0, 0
    ), nrow = 3, byrow = TRUE)
    expm::expm(Q * dt)
  }
  
  # helper: propagate one step from (H,U) using P
  step_forward <- function(H, U, P) {
    v_next <- c(H, U, 0) %*% P
    c(H = as.numeric(v_next[1]), U = as.numeric(v_next[2]))
  }
  
  # storage for hazards at each start-age
  hu <- rep(NA_real_, n)
  hd <- rep(NA_real_, n)
  ud <- rep(NA_real_, n)
  uh <- rep(min_haz, n)
  
  # solve interval-by-interval
  for (i in 1:(n - 1)) {
    
    H0 <- lh_tgt[i]
    U0 <- lu_tgt[i]
    
    # If everyone is dead (or almost), just carry tiny hazards
    if ((H0 + U0) <= min_haz) {
      hu[i] <- min_haz; hd[i] <- min_haz; ud[i] <- Rx[i] * hd[i]
      next
    }
    
    S1_tgt <- lx[i + 1]
    U1_tgt <- lu_tgt[i + 1]
    
    # inner solver: for given hu, find hd that matches S1_tgt
    solve_hd_given_hu <- function(hu_i) {
      f_hd <- function(hd_i) {
        P <- make_P(hu = hu_i, hd = hd_i, Rx_i = Rx[i], dt = age_int)
        nxt <- step_forward(H0, U0, P)
        (nxt["H"] + nxt["U"]) - S1_tgt
      }
      
      # bracket expansion if needed
      lo <- hd_bracket[1]; hi <- hd_bracket[2]
      f_lo <- f_hd(lo); f_hi <- f_hd(hi)
      
      k <- 0
      while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && k < 30) {
        hi <- hi * 2
        f_hi <- f_hd(hi)
        k <- k + 1
      }
      if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) return(NA_real_)
      
      uniroot(f_hd, lower = lo, upper = hi, tol = tol, maxiter = maxiter)$root
    }
    
    # outer solver: choose hu so that U1 matches U1_tgt (with hd chosen to match survival)
    f_hu <- function(hu_i) {
      hd_i <- solve_hd_given_hu(hu_i)
      if (!is.finite(hd_i)) return(NA_real_)
      
      P <- make_P(hu = hu_i, hd = hd_i, Rx_i = Rx[i], dt = age_int)
      nxt <- step_forward(H0, U0, P)
      nxt["U"] - U1_tgt
    }
    
    # bracket expansion for hu
    lo <- hu_bracket[1]; hi <- hu_bracket[2]
    f_lo <- f_hu(lo); f_hi <- f_hu(hi)
    
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
  
  # last age: 1-step extrapolation; these don't affect prev or lx matching.
  extrap_last <- function(x) {
    if (is.finite(x[n-1]) && is.finite(x[n-2]) && x[n-1] > 0 && x[n-2] > 0) {
      exp(2*log(x[n-1]) - log(x[n-2]))  # log-linear extrapolation
    } else if (is.finite(x[n-1])) {
      x[n-1]
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


# ---------
# (X) Final snap: smooth bump correction on hazards (typically HU/UH) to close HLE / integrand
# ---------

# Smooth bump basis on age midpoints; max=1 for stable coefficient scaling
bump_basis <- function(age, center = 75, width = 12, age_int = 1) {
  a_mid <- age + age_int/2
  g <- exp(-0.5 * ((a_mid - center) / width)^2)
  g / max(g)
}

# Apply multiplicative exp(c*g) correction to selected transitions
apply_bump_correction <- function(haz_long,
                                  coeffs,
                                  g,
                                  trans_allow = names(coeffs),
                                  clamp_coef = c(-1.5, 1.5),
                                  clamp_haz = c(1e-12, 5)) {
  
  if (is.null(names(coeffs))) stop("coeffs must be a named numeric vector.", call. = FALSE)
  
  # normalize hazards df
  haz_long <- haz_long |>
    dplyr::ungroup() |>
    dplyr::mutate(
      age    = as.integer(age),
      trans  = trimws(tolower(as.character(trans))),
      hazard = as.numeric(hazard)
    )
  
  # normalize coeffs + keep names
  nm <- trimws(tolower(names(coeffs)))
  coeffs <- pmin(clamp_coef[2], pmax(clamp_coef[1], coeffs))
  coeffs <- stats::setNames(as.numeric(coeffs), nm)
  # normalize allow list
  trans_allow <- trimws(tolower(trans_allow))
  
  ages <- sort(unique(haz_long$age))
  if (length(g) != length(ages)) stop("Length of g must equal number of unique ages in haz_long.", call. = FALSE)
  
  g_df <- tibble::tibble(age = ages, g = as.numeric(g))
  out <- dplyr::left_join(haz_long, g_df, by = "age")
  
  # explicit lookup: coefficient per row
  c_tr <- rep(0, nrow(out))
  ok <- out$trans %in% trans_allow
  if (any(ok)) {
    idx <- match(out$trans[ok], names(coeffs))
    c_tr[ok] <- coeffs[idx]
    c_tr[ok][is.na(c_tr[ok])] <- 0
  }
  
  bump <- c_tr * out$g
  out$hazard <- out$hazard * exp(bump)
  
  if (any(!is.finite(out$hazard))) {
    bad <- out[!is.finite(out$hazard), c("age","trans","hazard")][1:min(12, sum(!is.finite(out$hazard))), ]
    stop("apply_bump_correction produced non-finite hazards. Example rows:\n",
         paste(capture.output(print(bad)), collapse = "\n"),
         call. = FALSE)
  }
  
  out$hazard <- pmin(clamp_haz[2], pmax(clamp_haz[1], out$hazard))
  
  out |> dplyr::select(age, trans, hazard)
}
# Run MSLT directly from hazards (no need for theta)
run_mslt_from_hazards <- function(haz_long,
                                  age,
                                  init = c(H=1, U=0),
                                  age_int = 1) {
  
  age <- sort(age)
  
  # quick checks
  if (any(!is.finite(haz_long$hazard))) stop("Non-finite hazard values.", call. = FALSE)
  if (any(haz_long$hazard < 0)) stop("Negative hazards found.", call. = FALSE)
  
  lH <- numeric(length(age))
  lU <- numeric(length(age))
  lH[1] <- init[["H"]]
  lU[1] <- init[["U"]]
  
  hz_by_age <- split(haz_long, haz_long$age)
  
  for (i in seq_len(length(age) - 1)) {
    a <- as.character(age[i])
    hza <- hz_by_age[[a]]
    if (is.null(hza)) stop("Missing hazards at age ", a, call. = FALSE)
    
    hz_vec <- stats::setNames(hza$hazard, hza$trans)
    
    # require all 4 transitions
    need <- c("hu","hd","uh","ud")
    miss <- setdiff(need, names(hz_vec))
    if (length(miss) > 0) stop("Missing transitions at age ", a, ": ", paste(miss, collapse=", "), call. = FALSE)
    
    Q <- make_Q(hz_vec)
    
    # expm can error; let caller catch
    P <- expm::expm(Q * age_int)
    
    v <- c(H = lH[i], U = lU[i], D = 1 - (lH[i] + lU[i]))
    v_next <- as.numeric(v %*% P)
    names(v_next) <- c("H","U","D")
    
    lH[i+1] <- v_next[["H"]]
    lU[i+1] <- v_next[["U"]]
  }
  
  lx <- lH + lU
  prev <- ifelse(lx > 0, lU / lx, 0)
  
  tibble::tibble(age = age, lH = lH, lU = lU, lx = lx, prevalence = prev)
}
# Helper: Lx proxy from lx
lx_to_Lx <- function(lx, age_int = 1) {
  lx <- pmax(1e-12, lx)
  lx_next <- dplyr::lead(lx, default = dplyr::last(lx))
  (lx + lx_next) / 2 * age_int
}

# Objective for bump coefficients: match integrand and optionally prev/lx
bump_objective <- function(par,
                           haz_long_base,
                           g,
                           ground_summary,
                           age,
                           init = c(H=1,U=0),
                           age_int = 1,
                           trans = c("hu","uh"),
                           w_int = 10,
                           integrand = c("health","unhealthy"),
                           w_prev = 0,
                           w_lx = 0,
                           ridge = 1,
                           clamp_coef = c(-1.5, 1.5),
                           clamp_haz = c(1e-12, 5),
                           # ---- NEW ----
                           param_mode = c("free", "balance"),
                           delta_only = FALSE,
                           ridge_gamma = 10) {
  
  integrand <- match.arg(integrand)
  param_mode <- match.arg(param_mode)
  
  if (any(!is.finite(par))) return(1e12)
  
  # map -> per-transition coefficients
  coeffs <- bump_par_to_coeffs(par, trans = trans, param_mode = param_mode, delta_only = delta_only)
  
  # apply correction
  haz_adj <- tryCatch(
    apply_bump_correction(
      haz_long_base,
      coeffs = coeffs,
      g = g,
      trans_allow = trans,
      clamp_coef = clamp_coef,
      clamp_haz = clamp_haz
    ),
    error = function(e) NULL
  )
  if (is.null(haz_adj)) return(1e12)
  if (any(!is.finite(haz_adj$hazard))) return(1e12)
  
  # run MSLT; catch expm failures
  lt <- tryCatch(
    run_mslt_from_hazards(haz_adj, age = age, init = init, age_int = age_int),
    error = function(e) NULL
  )
  if (is.null(lt)) return(1e12)
  
  lt <- lt |>
    dplyr::left_join(
      ground_summary |> dplyr::select(age, lx_g = lx, prev_g = prevalence),
      by = "age"
    )
  
  # Lx proxies
  lx_to_Lx <- function(lx, age_int = 1) {
    lx <- pmax(1e-12, lx)
    lx_next <- dplyr::lead(lx, default = dplyr::last(lx))
    (lx + lx_next) / 2 * age_int
  }
  
  Lx1 <- lx_to_Lx(lt$lx, age_int)
  Lxg <- lx_to_Lx(lt$lx_g, age_int)
  
  if (integrand == "health") {
    I1 <- (1 - lt$prevalence) * Lx1
    Ig <- (1 - lt$prev_g) * Lxg
  } else {
    I1 <- lt$prevalence * Lx1
    Ig <- lt$prev_g * Lxg
  }
  
  loss_int  <- (I1 - Ig)^2
  loss_prev <- (lt$prevalence - lt$prev_g)^2
  loss_lx   <- (lt$lx - lt$lx_g)^2
  
  fit <- 0
  if (w_int  > 0) fit <- fit + w_int  * sum(loss_int,  na.rm = TRUE)
  if (w_prev > 0) fit <- fit + w_prev * sum(loss_prev, na.rm = TRUE)
  if (w_lx   > 0) fit <- fit + w_lx   * sum(loss_lx,   na.rm = TRUE)
  
  # ridge on parameters
  reg <- 0
  if (ridge > 0) reg <- reg + ridge * sum(as.numeric(par)^2)
  
  # extra penalty: keep gamma small (prevents “churn” solutions)
  if (param_mode == "balance" && !delta_only && ridge_gamma > 0) {
    gamma <- as.numeric(par[2])
    reg <- reg + ridge_gamma * gamma^2
  }
  
  fit + reg
}
# Solve bump coefficients for a given world (hazards)
calibrate_hazards_bump <- function(haz_long_base,
                                   ground_summary,
                                   age,
                                   init = c(H=1,U=0),
                                   age_int = 1,
                                   trans = c("hu","uh"),
                                   center = 75,
                                   width = 12,
                                   w_int = 10,
                                   integrand = c("health","unhealthy"),
                                   w_prev = 0,
                                   w_lx = 0,
                                   ridge = 1,
                                   clamp_coef = c(-1.5, 1.5),
                                   clamp_haz = c(1e-12, 5),
                                   method = c("optim","nlm"),
                                   maxit = 200,
                                   # ---- NEW ----
                                   param_mode = c("free", "balance"),
                                   delta_only = FALSE,
                                   ridge_gamma = 10) {
  
  integrand <- match.arg(integrand)
  method <- match.arg(method)
  param_mode <- match.arg(param_mode)
  
  g <- bump_basis(age, center = center, width = width, age_int = age_int)
  
  # parameter initialization
  if (param_mode == "free") {
    par0 <- rep(0, length(trans))
    names(par0) <- trans
  } else {
    # balance mode: delta (+ gamma unless delta_only)
    par0 <- if (delta_only) c(delta = 0) else c(delta = 0, gamma = 0)
  }
  
  fn <- function(p) bump_objective(
    par = p,
    haz_long_base = haz_long_base,
    g = g,
    ground_summary = ground_summary,
    age = age,
    init = init,
    age_int = age_int,
    trans = trans,
    w_int = w_int,
    integrand = integrand,
    w_prev = w_prev,
    w_lx = w_lx,
    ridge = ridge,
    clamp_coef = clamp_coef,
    clamp_haz = clamp_haz,
    param_mode = param_mode,
    delta_only = delta_only,
    ridge_gamma = ridge_gamma
  )
  
  if (method == "optim") {
    res <- stats::optim(
      par = par0,
      fn = fn,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-10)
    )
    par_hat <- res$par
    value   <- res$value
    conv    <- res$convergence
  } else {
    res <- stats::nlm(f = fn, p = par0, iterlim = maxit)
    par_hat <- res$estimate
    value   <- res$minimum
    conv    <- res$code
  }
  
  # map final params -> per-transition coefficients
  coeffs <- bump_par_to_coeffs(par_hat, trans = trans, param_mode = param_mode, delta_only = delta_only)
  
  haz_adj <- apply_bump_correction(
    haz_long_base,
    coeffs = coeffs,
    g = g,
    trans_allow = trans,
    clamp_coef = clamp_coef,
    clamp_haz = clamp_haz
  )
  
  list(
    par_hat = par_hat,     # <--- NEW: keep raw optimized params
    coeffs = coeffs,
    haz_adj = haz_adj,
    g = g,
    value = value,
    convergence = conv,
    param_mode = param_mode,
    delta_only = delta_only
  )
}
calibrate_hazards_multibump <- function(haz_long_base,
                                        ground_summary,
                                        age,
                                        init = c(H=1,U=0),
                                        age_int = 1,
                                        trans = c("hu","uh"),
                                        centers = seq(58, 80, by = 4),
                                        width = 6,
                                        w_int = 10,
                                        integrand = c("health","unhealthy"),
                                        w_prev = 0.5,
                                        w_lx = 0,
                                        ridge = 2,
                                        clamp_coef = c(-1, 1),
                                        clamp_haz = c(1e-12, 2),
                                        maxit = 200,
                                        method = "optim",
                                        # ---- NEW ----
                                        param_mode = c("free", "balance"),
                                        delta_only = FALSE,
                                        ridge_gamma = 10) {
  
  integrand <- match.arg(integrand)
  param_mode <- match.arg(param_mode)
  
  haz_cur <- haz_long_base
  path <- list()
  
  for (k in seq_along(centers)) {
    ck <- centers[k]
    
    cal_k <- calibrate_hazards_bump(
      haz_long_base = haz_cur,
      ground_summary = ground_summary,
      age = age,
      init = init,
      age_int = age_int,
      trans = trans,
      center = ck,
      width = width,
      w_int = w_int,
      integrand = integrand,
      w_prev = w_prev,
      w_lx = w_lx,
      ridge = ridge,
      clamp_coef = clamp_coef,
      clamp_haz = clamp_haz,
      method = method,
      maxit = maxit,
      param_mode = param_mode,
      delta_only = delta_only,
      ridge_gamma = ridge_gamma
    )
    
    haz_cur <- cal_k$haz_adj
    
    # record step info
    path[[k]] <- tibble::tibble(
      step = k,
      center = ck,
      width = width,
      param_mode = param_mode,
      delta_only = delta_only,
      value = cal_k$value,
      convergence = cal_k$convergence
    ) |>
      dplyr::bind_cols(as.data.frame(t(cal_k$par_hat))) |>
      dplyr::bind_cols(as.data.frame(t(cal_k$coeffs)))
  }
  
  list(
    haz_adj = haz_cur,
    path = dplyr::bind_rows(path)
  )
}
# Map bump optimization parameters to per-transition coefficients
# - free:   par is a vector of length length(trans), names = trans
# - balance: par is c(delta, gamma) or just c(delta) if delta_only=TRUE
bump_par_to_coeffs <- function(par,
                               trans,
                               param_mode = c("free", "balance"),
                               delta_only = FALSE) {
  param_mode <- match.arg(param_mode)
  
  trans <- trimws(tolower(trans))
  
  if (param_mode == "free") {
    coeffs <- stats::setNames(as.numeric(par), trans)
    return(coeffs)
  }
  
  # balance mode requires exactly two transitions
  if (length(trans) != 2) stop("param_mode='balance' requires exactly two transitions in `trans`.", call. = FALSE)
  
  if (delta_only) {
    if (length(par) != 1) stop("delta_only=TRUE expects par length 1 (delta).", call. = FALSE)
    delta <- as.numeric(par[1])
    gamma <- 0
  } else {
    if (length(par) != 2) stop("param_mode='balance' expects par length 2 (delta, gamma).", call. = FALSE)
    delta <- as.numeric(par[1])
    gamma <- as.numeric(par[2])
  }
  
  c1 <-  delta + gamma
  c2 <- -delta + gamma
  
  stats::setNames(c(c1, c2), trans)
}



# --- helpers to get targets from ground_summary ---
ground_to_lH_lU <- function(ground_summary) {
  gs <- ground_summary
  if (!all(c("age","lx","prevalence") %in% names(gs))) {
    stop("ground_summary must have columns: age, lx, prevalence", call. = FALSE)
  }
  gs <- gs |>
    dplyr::mutate(
      age = as.integer(age),
      lH  = lx * (1 - prevalence),
      lU  = lx * prevalence
    )
  gs
}

# --- one-step projection given hazards at age x ---
step_state <- function(lH, lU, haz_vec, age_int = 1) {
  Q <- make_Q(haz_vec)
  P <- expm::expm(Q * age_int)
  v <- c(H = lH, U = lU, D = 1 - (lH + lU))
  v_next <- as.numeric(v %*% P)
  c(H = v_next[1], U = v_next[2], D = v_next[3])
}

# --- objective for a single age ---
snap_obj_one_age <- function(par_log,
                             trans,                 # transitions being optimized
                             lH, lU,
                             targetH, targetU,
                             haz_base,              # named hu,hd,uh,ud baseline hazards at age
                             age_int = 1,
                             wH = 1, wU = 1,
                             # optional continuity penalty (on log-scale)
                             prev_par_log = NULL,   # named vector for same trans from previous age (log-haz)
                             smooth_w = 0) {
  
  if (any(!is.finite(par_log))) return(1e12)
  
  # build hazards for this age
  haz <- haz_base
  haz[trans] <- exp(as.numeric(par_log))
  
  pred <- tryCatch(step_state(lH, lU, haz, age_int), error = function(e) NULL)
  if (is.null(pred) || any(!is.finite(pred))) return(1e12)
  
  loss <- wH * (pred[["H"]] - targetH)^2 + wU * (pred[["U"]] - targetU)^2
  
  # smoothness penalty on log hazards (preferred)
  if (smooth_w > 0 && !is.null(prev_par_log)) {
    prev_par_log <- prev_par_log[trans]
    ok <- is.finite(prev_par_log)
    if (any(ok)) {
      loss <- loss + smooth_w * sum((as.numeric(par_log)[ok] - as.numeric(prev_par_log)[ok])^2)
    }
  }
  
  loss
}
# --- main forward loop: calibrate hazards age-by-age ---
snap_hazards_agewise <- function(haz_long_base,
                                 ground_summary,
                                 age,
                                 age_int = 1,
                                 trans = c("hu"),
                                 # bounds are on hazards, but converted to log-bounds
                                 bounds = NULL,            # NULL / c(lo,hi) / named list per trans
                                 wH = 1, wU = 1,
                                 smooth_w = 0,
                                 maxit = 400,
                                 method = c("auto", "L-BFGS-B", "BFGS", "Nelder-Mead"),
                                 verbose = FALSE) {
  
  method <- match.arg(method)
  age <- sort(as.integer(age))
  trans <- trimws(tolower(trans))
  need_all <- c("hu","hd","uh","ud")
  if (!all(trans %in% need_all)) stop("trans must be subset of c('hu','hd','uh','ud')", call. = FALSE)
  if (length(trans) < 1 || length(trans) > 4) stop("trans must have length 1..4", call. = FALSE)
  
  # normalize hazards input
  haz0 <- haz_long_base |>
    dplyr::ungroup() |>
    dplyr::mutate(
      age = as.integer(age),
      trans = trimws(tolower(as.character(trans))),
      hazard = as.numeric(hazard)
    ) |>
    dplyr::arrange(age, trans)
  
  hz_by_age <- split(haz0, haz0$age)
  
  # ground targets
  gs <- ground_to_lH_lU(ground_summary)
  lH_ground <- stats::setNames(gs$lH, gs$age)
  lU_ground <- stats::setNames(gs$lU, gs$age)
  
  # bounds handling (in hazard space)
  if (is.null(bounds)) {
    bounds <- setNames(rep(list(c(1e-12, 2)), length(trans)), trans)
  } else {
    if (is.numeric(bounds) && length(bounds) == 2 && is.null(names(bounds))) {
      bounds <- setNames(rep(list(bounds), length(trans)), trans)
    } else if (is.list(bounds)) {
      if (is.null(names(bounds))) stop("bounds list must be named by transition.", call. = FALSE)
      for (tr in trans) if (is.null(bounds[[tr]])) bounds[[tr]] <- c(1e-12, 2)
      bounds <- bounds[trans]
    } else {
      stop("bounds must be NULL, a length-2 numeric vector, or a named list.", call. = FALSE)
    }
  }
  
  lower_h <- vapply(trans, function(tr) bounds[[tr]][1], numeric(1))
  upper_h <- vapply(trans, function(tr) bounds[[tr]][2], numeric(1))
  
  # convert to log-bounds for L-BFGS-B
  lower_log <- log(pmax(1e-300, lower_h))
  upper_log <- log(pmax(1e-300, upper_h))
  
  # output hazards updated age-by-age
  haz_out <- haz0
  
  # continuity penalty uses prev log-params
  prev_par_log <- NULL
  
  for (i in seq_len(length(age) - 1)) {
    a  <- age[i]
    a1 <- age[i+1]
    
    hza <- hz_by_age[[as.character(a)]]
    if (is.null(hza)) stop("Missing hazards at age ", a, call. = FALSE)
    
    haz_base <- stats::setNames(hza$hazard, hza$trans)
    miss <- setdiff(need_all, names(haz_base))
    if (length(miss) > 0) stop("Missing transitions at age ", a, ": ", paste(miss, collapse=", "), call. = FALSE)
    
    # start/target (strictly from ground world per your spec)
    lH <- unname(lH_ground[[as.character(a)]])
    lU <- unname(lU_ground[[as.character(a)]])
    targetH <- unname(lH_ground[[as.character(a1)]])
    targetU <- unname(lU_ground[[as.character(a1)]])
    
    if (!all(is.finite(c(lH,lU,targetH,targetU)))) stop("Non-finite ground lH/lU at age ", a, call. = FALSE)
    
    # starting point in log space
    par_start_log <- log(pmax(1e-300, as.numeric(haz_base[trans])))
    names(par_start_log) <- trans
    
    fn <- function(p_log) snap_obj_one_age(
      par_log = p_log,
      trans = trans,
      lH = lH, lU = lU,
      targetH = targetH, targetU = targetU,
      haz_base = haz_base,
      age_int = age_int,
      wH = wH, wU = wU,
      prev_par_log = prev_par_log,
      smooth_w = smooth_w
    )
    
    # choose method automatically if requested
    meth <- method
    if (meth == "auto") {
      meth <- if (length(trans) == 1) "optimize" else "L-BFGS-B"
    }
    
    if (meth == "optimize") {
      if (length(trans) != 1) stop("method='optimize' only works when length(trans)==1", call. = FALSE)
      f1 <- function(x) fn(x)
      res <- stats::optimize(f1, interval = c(lower_log[1], upper_log[1]))
      par_hat_log <- res$minimum
      names(par_hat_log) <- trans
      
    } else if (meth == "L-BFGS-B") {
      res <- stats::optim(
        par = par_start_log,
        fn  = fn,
        method = "L-BFGS-B",
        lower = lower_log,
        upper = upper_log,
        control = list(maxit = maxit)
      )
      par_hat_log <- res$par
      names(par_hat_log) <- trans
      
    } else {
      # unconstrained: BFGS / Nelder-Mead
      # NOTE: without bounds, hazards can get extreme; clamp via a soft penalty:
      fn_uncon <- function(p_log) {
        # soft penalty if outside log-bounds
        pen <- 0
        if (any(p_log < lower_log) || any(p_log > upper_log)) {
          d <- pmax(0, lower_log - p_log) + pmax(0, p_log - upper_log)
          pen <- 1e6 * sum(d^2)
        }
        fn(p_log) + pen
      }
      
      res <- stats::optim(
        par = par_start_log,
        fn  = fn_uncon,
        method = meth,
        control = list(maxit = maxit)
      )
      par_hat_log <- res$par
      names(par_hat_log) <- trans
      # clamp result back into bounds
      par_hat_log <- pmin(upper_log, pmax(lower_log, par_hat_log))
    }
    
    # update hazards at age a
    haz_base[trans] <- exp(as.numeric(par_hat_log))
    
    idx <- haz_out$age == a
    haz_out$hazard[idx] <- haz_base[haz_out$trans[idx]]
    
    # update previous log-params for smoothness penalty
    prev_par_log <- stats::setNames(as.numeric(par_hat_log), trans)
    
    if (verbose && (i %% 5 == 0)) {
      pred <- step_state(lH, lU, haz_base, age_int)
      cat("age", a, "loss", (pred["H"]-targetH)^2 + (pred["U"]-targetU)^2, "\n")
    }
  }
  
  haz_out
}


#----------------------------------------------------
# ---- numerical Jacobian (finite differences on log-scale for stability) ----
jac_fd_log <- function(f, z, eps = 1e-5) {
  # z is log-parameters
  f0 <- f(z)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(z))
  for (j in seq_along(z)) {
    zp <- z; zp[j] <- zp[j] + eps
    zm <- z; zm[j] <- zm[j] - eps
    fp <- f(zp)
    fm <- f(zm)
    J[, j] <- (fp - fm) / (2 * eps)
  }
  J
}

# ---- damped Newton in log-parameter space ----
newton_root_log2 <- function(f, z0,
                             tol = 1e-10,
                             maxit = 50,
                             jac_eps = 1e-5,
                             step_max = 1.0,
                             verbose = FALSE) {
  z <- z0
  r <- f(z)
  val <- sum(r^2)
  
  for (it in seq_len(maxit)) {
    if (!all(is.finite(r))) return(list(converged = FALSE, z = z, it = it, val = Inf))
    
    if (sqrt(val) < tol) return(list(converged = TRUE, z = z, it = it, val = val))
    
    J <- jac_fd_log(f, z, eps = jac_eps)
    
    # Solve J * dz = -r (least squares if singular)
    dz <- tryCatch(
      as.numeric(qr.solve(J, -r)),
      error = function(e) {
        # fallback: ridge-stabilized normal equations
        JTJ <- t(J) %*% J
        as.numeric(solve(JTJ + diag(1e-8, ncol(JTJ)), t(J) %*% (-r)))
      }
    )
    
    # clamp step size (helps avoid wild swings)
    dz_norm <- sqrt(sum(dz^2))
    if (dz_norm > step_max) dz <- dz * (step_max / dz_norm)
    
    # backtracking line search
    tstep <- 1
    ok <- FALSE
    for (ls in 1:20) {
      z_try <- z + tstep * dz
      r_try <- f(z_try)
      val_try <- sum(r_try^2)
      
      if (is.finite(val_try) && val_try < val) {
        z <- z_try
        r <- r_try
        val <- val_try
        ok <- TRUE
        break
      }
      tstep <- tstep / 2
    }
    
    if (!ok) {
      if (verbose) cat("Newton failed line search at it=", it, "val=", val, "\n")
      return(list(converged = FALSE, z = z, it = it, val = val))
    }
    
    if (verbose) cat("it", it, "val", signif(val, 6), "\n")
  }
  
  list(converged = FALSE, z = z, it = maxit, val = val)
}

# ---- one-age root solve for (hu, uh) with fixed (hd, ud) ----
solve_hu_uh_one_age <- function(lH, lU, targetH, targetU,
                                hd, ud,
                                hu0, uh0,
                                age_int = 1,
                                lower = c(1e-12, 1e-12),
                                upper = c(5, 5),
                                tol = 1e-10,
                                maxit = 50,
                                verbose = FALSE) {
  
  # work on log scale with bounds enforced by clamping in f()
  lo <- log(lower); hi <- log(upper)
  z0 <- log(pmax(lower, pmin(upper, c(hu0, uh0))))
  
  f <- function(z) {
    z <- pmin(hi, pmax(lo, z))
    hu <- exp(z[1]); uh <- exp(z[2])
    haz <- c(hu = hu, hd = hd, uh = uh, ud = ud)
    pred <- step_state(lH, lU, haz, age_int = age_int)
    c(pred[["H"]] - targetH,
      pred[["U"]] - targetU)
  }
  
  res <- newton_root_log2(f, z0, tol = tol, maxit = maxit, verbose = verbose)
  
  zhat <- pmin(hi, pmax(lo, res$z))
  c(hu = exp(zhat[1]), uh = exp(zhat[2]), converged = res$converged, val = res$val)
}




