# State-specific lifespan variation across Sullivan-compatible worlds
#
# Outputs:
#   figs/si_variance_components_monthly.pdf
#   figs/si_sullivan_validation_monthly.pdf
#   tables/si_variance_ranges_monthly.tex
#
# Change age_int to 1 for the corresponding annual sensitivity outputs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(xtable)
})

# -----------------------------------------------------------------------------
# Settings and data
# -----------------------------------------------------------------------------

age_int <- 1 / 12
init <- c(H = 1, U = 0)

interval_label <- dplyr::case_when(
  isTRUE(all.equal(age_int, 1))      ~ "annual",
  isTRUE(all.equal(age_int, 1 / 12)) ~ "monthly",
  TRUE                               ~ paste0("step_", format(age_int))
)

data_file <- dplyr::case_when(
  interval_label == "annual"  ~ "data/worlds_mslt_annual.csv.gz",
  interval_label == "monthly" ~ "data/worlds_mslt_monthly.csv.gz",
  TRUE ~ NA_character_
)

if (is.na(data_file)) {
  stop("No data file has been specified for age_int = ", age_int, call. = FALSE)
}

dir.create("figs", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

worlds <- readr::read_csv(data_file, show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Multistate joint distribution and moments
# -----------------------------------------------------------------------------

# Efficient recursion over age and accumulated healthy duration. This follows
# the discrete-time accounting in Riffe et al. (2024): at exact duration x,
# accumulated durations satisfy x = h + u. The terminal age is closed by
# assigning all remaining survivors to death.
calc_dxh_fast <- function(p_tibble,
                          age_int = 1,
                          init = c(H = 1, U = 0),
                          keep_zero = FALSE) {
  
  p_tibble <- p_tibble |>
    dplyr::arrange(age)
  
  ages <- p_tibble$age
  min_age <- min(ages)
  n_ages <- length(ages)
  
  H <- unname(init["H"])
  U <- unname(init["U"])
  out <- vector("list", n_ages)
  
  for (i in seq_len(n_ages)) {
    age_i <- ages[i]
    x_i <- age_i - min_age
    h_i <- seq.int(0, length(H) - 1L) * age_int
    u_i <- x_i - h_i
    terminal <- i == n_ages
    
    qH <- if (terminal) 1 else p_tibble$HD[i]
    qU <- if (terminal) 1 else p_tibble$UD[i]
    
    out[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        current_state = "H",
        age = age_i,
        x = x_i,
        h = h_i,
        u = u_i,
        lxsc = H,
        dxsc = H * qH
      ),
      tibble::tibble(
        current_state = "U",
        age = age_i,
        x = x_i,
        h = h_i,
        u = u_i,
        lxsc = U,
        dxsc = U * qU
      )
    )
    
    if (!terminal) {
      # Arrival in H increments accumulated healthy duration by one step.
      H_next <- numeric(length(H) + 1L)
      H_next[-1L] <- H * p_tibble$HH[i] + U * p_tibble$UH[i]
      
      # Arrival in U leaves accumulated healthy duration unchanged.
      U_next <- numeric(length(U) + 1L)
      U_next[seq_along(U)] <- U * p_tibble$UU[i] + H * p_tibble$HU[i]
      
      H <- H_next
      U <- U_next
    }
  }
  
  ans <- dplyr::bind_rows(out)
  if (!keep_zero) {
    ans <- dplyr::filter(ans, lxsc > 0 | dxsc > 0)
  }
  ans
}

joint_moments <- function(p_tibble,
                          age_int = 1,
                          init = c(H = 1, U = 0)) {
  
  d <- calc_dxh_fast(
    p_tibble = p_tibble,
    age_int = age_int,
    init = init
  ) |>
    dplyr::filter(dxsc > 0) |>
    dplyr::mutate(probability = dxsc / sum(dxsc))
  
  # Do not add a half-interval correction here. The multistate distribution is
  # indexed on its discrete duration coordinates, for which x = h + u exactly.
  le  <- sum(d$x * d$probability)
  hle <- sum(d$h * d$probability)
  ule <- sum(d$u * d$probability)
  
  vle <- sum((d$x - le)^2 * d$probability)
  vh  <- sum((d$h - hle)^2 * d$probability)
  vu  <- sum((d$u - ule)^2 * d$probability)
  cov_hu <- sum((d$h - hle) * (d$u - ule) * d$probability)
  
  tibble::tibble(
    le = le,
    hle = hle,
    ule = ule,
    vle = vle,
    vh = vh,
    vu = vu,
    cov_hu = cov_hu,
    vle_check = vh + vu + 2 * cov_hu,
    check_difference = vle - vle_check,
    probability_check = sum(d$probability)
  )
}

# Summarize within each group immediately. This avoids retaining the full joint
# distribution for all 32 worlds in memory.
variances <- worlds |>
  dplyr::select(system, world, age, HH, HU, HD, UH, UU, UD) |>
  dplyr::group_by(system, world) |>
  dplyr::group_modify(
    ~ joint_moments(.x, age_int = age_int, init = init)
  ) |>
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# Sullivan healthy-lifespan distribution
# -----------------------------------------------------------------------------

# Muszynska-Spielauer et al. (2024), especially equations (9), (12),
# and (15)-(16). Rows represent age intervals [x, x+n).
sullivan_healthy_dist <- function(age,
                                  Lx,
                                  prevalence,
                                  lx = NULL,
                                  n = NULL,
                                  ax = NULL,
                                  terminal = c("close", "drop"),
                                  tol = 1e-10) {
  
  terminal <- match.arg(terminal)
  N <- length(age)
  stopifnot(length(Lx) == N, length(prevalence) == N)
  
  if (is.null(n)) {
    d_age <- diff(age)
    if (!length(d_age)) stop("At least two ages are required.")
    n <- c(d_age, tail(d_age, 1L))
  }
  if (length(n) == 1L) n <- rep(n, N)
  if (is.null(ax)) ax <- n / 2
  if (length(ax) == 1L) ax <- rep(ax, N)
  stopifnot(length(n) == N, length(ax) == N, all(n > 0))
  
  # Healthy stationary-population stock, interpreted as healthy survivorship
  # at exact age x + ax under stationarity and irreversibility.
  healthy_stock <- (1 - prevalence) * Lx / n
  increase <- c(FALSE, diff(healthy_stock) > tol)
  
  if (any(increase)) {
    stop(
      "Healthy stock increases with age; an irreversible healthy-lifespan ",
      "distribution is not identified. First offending age: ",
      age[which(increase)[1L]],
      call. = FALSE
    )
  }
  
  exits <- pmax(
    healthy_stock - dplyr::lead(healthy_stock, default = 0),
    0
  )
  if (terminal == "drop") exits[N] <- 0
  
  age_ax <- age + ax
  exit_age <- (
    age_ax + dplyr::lead(age_ax, default = age_ax[N])
  ) / 2
  
  if (!is.null(lx)) {
    stopifnot(length(lx) == N)
    dx <- pmax(lx - dplyr::lead(lx, default = 0), 0)
    d_second <- dx * ax / n
    d_first_next <- dplyr::lead(dx * (1 - ax / n), default = 0)
    den <- d_second + d_first_next
    use <- seq_len(N - 1L)
    
    exit_age[use] <- ifelse(
      den[use] > 0,
      (
        d_second[use] * age_ax[use] +
          d_first_next[use] * age_ax[use + 1L]
      ) / den[use],
      exit_age[use]
    )
    exit_age[N] <- age_ax[N]
  }
  
  radix <- healthy_stock[1L]
  
  tibble::tibble(
    age = age,
    age_ax = age_ax,
    healthy_stock = healthy_stock,
    exit_age = exit_age,
    healthy_duration = exit_age - age[1L],
    exits = exits,
    probability = exits / radix
  )
}

dist_moments <- function(duration, probability) {
  keep <- is.finite(duration) &
    is.finite(probability) &
    probability > 0
  
  duration <- duration[keep]
  probability <- probability[keep]
  probability <- probability / sum(probability)
  
  mean_duration <- sum(duration * probability)
  variance <- sum((duration - mean_duration)^2 * probability)
  
  tibble::tibble(
    mean = mean_duration,
    variance = variance,
    sd = sqrt(variance)
  )
}

compare_sullivan_variance <- function(worlds, age_int) {
  worlds |>
    dplyr::group_by(system, world) |>
    dplyr::group_modify(~ {
      x <- dplyr::arrange(.x, age)
      d <- sullivan_healthy_dist(
        age = x$age,
        Lx = x$Lx,
        prevalence = x$prevalence_interval,
        lx = x$lx,
        n = age_int,
        ax = age_int / 2
      )
      dist_moments(d$healthy_duration, d$probability)
    }) |>
    dplyr::ungroup()
}

sullivan_variance <- compare_sullivan_variance(worlds, age_int)

# It is mechanically common to all worlds because their Sullivan inputs are
# common, but it has a healthy-lifespan interpretation only for no-returns worlds.
sullivan_no_returns <- sullivan_variance |>
  dplyr::filter(system == "noreturns")

validation <- variances |>
  dplyr::filter(system == "noreturns") |>
  dplyr::select(system, world, multistate_variance = vh) |>
  dplyr::left_join(
    dplyr::select(
      sullivan_no_returns,
      system,
      world,
      sullivan_variance = variance
    ),
    by = c("system", "world")
  ) |>
  dplyr::mutate(
    difference = sullivan_variance - multistate_variance
  )

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

diagnostics <- variances |>
  dplyr::summarise(
    max_probability_error = max(abs(probability_check - 1)),
    max_variance_identity_error = max(abs(check_difference)),
    .by = system
  )

no_returns_checks <- variances |>
  dplyr::filter(system == "noreturns") |>
  dplyr::summarise(
    healthy_variance_range = max(vh) - min(vh),
    unhealthy_variance_range = max(vu) - min(vu),
    covariance_range = max(cov_hu) - min(cov_hu),
    total_variance_range = max(vle) - min(vle)
  )

print(diagnostics)
print(no_returns_checks)
print(validation)

# -----------------------------------------------------------------------------
# Supplementary Figure S5: state-specific variance components
# -----------------------------------------------------------------------------

system_labels <- c(
  returns   = "Returns",
  noreturns = "No returns"
)

measure_labels <- c(
  vh     = "Var(H)",
  vu     = "Var(U)",
  cov_hu = "Cov(H,U)"
)

measure_colors <- c(
  "Var(H)"   = "#00BFC4",
  "Var(U)"   = "#F8766D",
  "Cov(H,U)" = "#C77CFF"
)

component_data <- variances |>
  select(system, world, vh, vu, cov_hu) |>
  pivot_longer(
    cols = c(vh, vu, cov_hu),
    names_to = "measure",
    values_to = "value"
  ) |>
  mutate(
    system = factor(
      system,
      levels = names(system_labels),
      labels = unname(system_labels)
    ),
    measure = factor(
      measure,
      levels = names(measure_labels),
      labels = unname(measure_labels)
    ),
    # World 1 at the top and world 16 at the bottom.
    world = factor(
      world,
      levels = rev(sort(unique(world)))
    )
  )

p_s5 <- ggplot(
  component_data,
  aes(
    x = value,
    y = world,
    color = measure
  )
) +
  geom_point(
    size = 2.8,
    alpha = 0.9
  ) +
  facet_grid(
    rows = vars(system),
    cols = vars(measure),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = measure_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    x = expression(paste("Years"^2)),
    y = "World"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.35
    ),
    panel.grid.major.x = element_blank(),
    
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 14
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 14,
      angle = 0
    ),
    
    axis.text = element_text(size = 10.5),
    axis.title = element_text(size = 13),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      margin = margin(r = 10)
    ),
    
    panel.spacing.x = unit(1.5, "lines"),
    panel.spacing.y = unit(1.1, "lines"),
    plot.margin = margin(8, 12, 8, 8)
  )

fig_s5_file <- file.path(
  "figs",
  paste0(
    "si_variance_components_",
    interval_label,
    ".pdf"
  )
)

ggsave(
  filename = fig_s5_file,
  plot = p_s5,
  width = 24,
  height = 15,
  units = "cm",
  device = cairo_pdf
)


# -----------------------------------------------------------------------------
# Supplementary Table S1: ranges across worlds
# -----------------------------------------------------------------------------

table_s1 <- variances |>
  dplyr::select(system, vh, vu, cov_hu, vle) |>
  tidyr::pivot_longer(
    cols = c(vh, vu, cov_hu, vle),
    names_to = "quantity",
    values_to = "value"
  ) |>
  dplyr::summarise(
    Minimum = min(value),
    Maximum = max(value),
    .by = c(system, quantity)
  ) |>
  dplyr::mutate(
    System = dplyr::recode(
      system,
      returns = "Returns",
      noreturns = "No returns"
    ),
    Quantity = dplyr::recode(
      quantity,
      vh = "$\\operatorname{Var}(H)$",
      vu = "$\\operatorname{Var}(U)$",
      cov_hu = "$\\operatorname{Cov}(H,U)$",
      vle = "$\\operatorname{Var}(X)$"
    ),
    quantity_order = match(quantity, c("vh", "vu", "cov_hu", "vle")),
    system_order = match(system, c("noreturns", "returns"))
  ) |>
  dplyr::arrange(system_order, quantity_order) |>
  dplyr::select(System, Quantity, Minimum, Maximum)

sullivan_table_row <- tibble::tibble(
  System = "No returns (Sullivan)",
  Quantity = "$\\operatorname{Var}(H)$",
  Minimum = min(sullivan_no_returns$variance),
  Maximum = max(sullivan_no_returns$variance)
)

table_s1 <- dplyr::bind_rows(table_s1, sullivan_table_row)

table_s1_xtable <- xtable::xtable(
  table_s1,
  caption = paste(
    "Ranges of state-specific lifespan variances and covariance across",
    "the 16 worlds in each transition system. The Sullivan row gives the",
    "healthy-years variance obtained from the common Sullivan inputs using",
    "the no-recovery approximation."
  ),
  label = "tab:variance_ranges",
  align = c("l", "l", "l", "r", "r"),
  digits = c(0, 0, 0, 3, 3)
)

table_s1_file <- file.path(
  "tables",
  paste0("si_variance_ranges_", interval_label, ".tex")
)

print(
  table_s1_xtable,
  file = table_s1_file,
  include.rownames = FALSE,
  floating = TRUE,
  table.placement = "htbp",
  caption.placement = "top",
  sanitize.text.function = identity,
  hline.after = c(-1, 0, nrow(table_s1))
)


