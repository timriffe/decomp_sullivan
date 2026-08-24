age_int <- 1 / 12
library(tidyverse)

ground_mslt <- read_csv("data/ground_mslt_monthly.csv")

str(ground_mslt)
# -----------------------------------------------------------------------------
# Supplementary Figure S1: ground-world Sullivan inputs
# -----------------------------------------------------------------------------

fig_s1_data <- ground_mslt |>
  transmute(age, Survivorship = lx, `Point prevalence` = prevalence_point) |>
  pivot_longer(cols = -age,
               names_to = "measure",
               values_to = "value") |>
  mutate(measure = factor(measure, levels = c("Survivorship", "Point prevalence")))

fig_s1_labels <- c("Survivorship  \u2113(a)", "Prevalence  \u03C0(a)")

p_s1 <- ggplot(fig_s1_data, aes(x = age, y = value, color = measure)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Survivorship" = "#00BFC4",
      "Point prevalence" = "#F8766D"
    ),
    breaks = c("Survivorship", "Point prevalence"),
    labels = fig_s1_labels
  ) +
  scale_x_continuous(breaks = seq(50, 100, by = 10), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Age", y = "Proportion", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.title.x = element_text(margin = margin(t = 9)),
    axis.title.y = element_text(margin = margin(r = 9)),
    
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.width = unit(1.8, "lines"),
    
    plot.margin = margin(8, 10, 8, 8)
  )

fig_s1_file <- file.path("figs", "si_ground_world_monthly.pdf")

ggsave(
  filename = fig_s1_file,
  plot = p_s1,
  width = 16,
  height = 10.5,
  units = "cm",
  device = cairo_pdf
)

p_s1

# Figure S2: all hazards:

haz <- read_csv("data/worlds_hazards_monthly.csv.gz")

# -----------------------------------------------------------------------------
# Supplementary Figure S2: hazards in all 32 simulated worlds
# -----------------------------------------------------------------------------

transition_levels <- c("hu", "uh", "hd", "ud")

transition_labels <- c(hu = "h[hu]",
                       uh = "h[uh]",
                       hd = "h[hd]",
                       ud = "h[ud]")

system_colors <- c("Returns"    = "#00BFC4",
                   "No returns" = "#F8766D")

haz_s2 <- haz |>
  filter(system %in% c("returns", "noreturns"), world %in% 1:16) |>
  mutate(
    system = factor(
      system,
      levels = c("returns", "noreturns"),
      labels = c("Returns", "No returns")
    ),
    trans = factor(
      trans,
      levels = transition_levels,
      labels = unname(transition_labels)
    ),
    world = factor(world, levels = 1:16)
  )

# Remove the numerical floor used to represent zero recovery.
haz_s2_lines <- haz_s2 |>
  filter(!(system == "No returns" &
             trans == "h[uh]"))

# Put the zero-recovery annotation near the middle of the vertical range
# occupied by the actual returns-system recovery hazards.
uh_label_y <- haz |>
  filter(system == "returns", trans == "uh", hazard > 0) |>
  summarize(value = exp(mean(log(hazard)))) |>
  pull(value)

zero_recovery_label <- tibble(
  age = 75,
  hazard = uh_label_y,
  system = factor("No returns", levels = c("Returns", "No returns")),
  trans = factor("h[uh]", levels = unname(transition_labels)),
  label = "0 by construction"
)

p_s2 <- ggplot(haz_s2_lines, aes(
  x = age,
  y = hazard,
  group = world,
  color = system
)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_text(
    data = zero_recovery_label,
    aes(x = age, y = hazard, label = label),
    inherit.aes = FALSE,
    color = "grey35",
    size = 4
  ) +
  facet_grid(
    rows = vars(trans),
    cols = vars(system),
    scales = "free_y",
    drop = FALSE,
    switch = "y",
    labeller = labeller(trans = label_parsed)
  ) +
  scale_color_manual(values = system_colors, guide = "none") +
  scale_x_continuous(breaks = seq(50, 100, by = 10), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_log10(
    labels = function(x) {
      sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
    }
  ) +
  labs(x = "Age", y = "Hazard") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    
    
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.0, "lines"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y.left = element_text(angle = 0, size = 13),
    
    axis.text = element_text(size = 10.5),
    axis.text.y = element_text(hjust = 0),
    axis.title = element_text(size = 13),
    axis.title.x = element_text(margin = margin(t = 9)),
    axis.title.y = element_text(margin = margin(r = 9)),
    
    plot.margin = margin(8, 10, 8, 8)
  )

fig_s2_file <- file.path("figs", "si_all_world_hazards_monthly.pdf")

ggsave(
  filename = fig_s2_file,
  plot = p_s2,
  width = 20,
  height = 22,
  units = "cm",
  device = cairo_pdf
)

p_s2


# =============================================================================
# Supplementary Figure S3:
# transition-specific decomposition margins, all 16 worlds
# =============================================================================

# Expects:
#   age_int <- 1       or
#   age_int <- 1 / 12
#

interval_label <- if_else(age_int == 1, "annual", "monthly")

# -----------------------------------------------------------------------------
# Read and prepare decomposition results
# -----------------------------------------------------------------------------

if (age_int == 1) {
  dec_s3 <- read_csv("data/worlds_decomp_annual.csv.gz", show_col_types = FALSE)
}

if (age_int == 1 / 12) {
  dec_s3 <- read_csv("data/worlds_decomp_monthly.csv.gz", show_col_types = FALSE) |>
    # Aggregate monthly contributions to completed single ages.
    mutate(age = age - age %% 1) |>
    summarize(cc = sum(cc, na.rm = TRUE),
              .by = c(system, world1, world2, trans, age))
}

# Preserve the order and colors in the existing acceptable figures:
# hd = coral, hu = green, ud = cyan, uh = purple.
transition_levels_s3 <- c("hd", "hu", "ud", "uh")

transition_colors_s3 <- c(
  hd = "#F8766D",
  hu = "#7CAE00",
  ud = "#00BFC4",
  uh = "#C77CFF"
)

transition_labels_s3 <- parse(text = c("h[hd]", "h[hu]", "h[ud]", "h[uh]"))

# ------------------------------------------------------------
# Add the reverse of every observed comparison
# ------------------------------------------------------------
# ------------------------------------------------------------
# Add the reverse of every observed comparison
# ------------------------------------------------------------

dec_s3_forward <-
  dec_s3 |>
  mutate(
    world1 = as.integer(as.character(world1)),
    world2 = as.integer(as.character(world2))
  )

dec_s3_reverse <-
  dec_s3_forward |>
  rename(
    world1_original = world1,
    world2_original = world2
  ) |>
  transmute(
    system,
    world1 = world2_original,
    world2 = world1_original,
    trans,
    age,
    cc = -cc
  )

dec_s3_symmetric <-
  bind_rows(dec_s3_forward, dec_s3_reverse) |>
  filter(world1 != world2) |>
  mutate(
    world1 = factor(world1, levels = 1:16),
    world2 = factor(world2, levels = 1:16)
  )

# Sum age-specific contributions into transition margins
dec_s3_margins <-
  dec_s3_symmetric |>
  group_by(system, world1, world2, trans) |>
  summarise(
    cc = sum(cc, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Function for one transition system
# -----------------------------------------------------------------------------

make_s3_margin_plot <- function(system_value) {
  plot_data <- dec_s3_margins |>
    filter(system == system_value)
  
  # Round upward to the next half-year. This gives approximately ±1.5
  # for returns and ±3 for no returns with the current simulations.
  plot_limit <- max(abs(plot_data$cc), na.rm = TRUE)
  
  plot_limit <- ceiling(plot_limit * 2) / 2
  
  if (!is.finite(plot_limit) ||
      plot_limit == 0) {
    plot_limit <- 0.5
  }
  
  # All off-diagonal comparison panels plus the upper-left panel,
  # which is reserved for the small reference axis.
  panel_grid <- expand_grid(world1 = factor(1:16, levels = 1:16),
                            world2 = factor(1:16, levels = 1:16)) |>
    filter(world1 != world2 |
             (
               world1 == factor(1, levels = 1:16) &
                 world2 == factor(1, levels = 1:16)
             ))
  
  background_negative <- panel_grid |>
    mutate(
      xmin = 0,
      xmax = 4.5,
      ymin = -Inf,
      ymax = 0
    )
  
  background_positive <- panel_grid |>
    mutate(
      xmin = 0,
      xmax = 4.5,
      ymin = 0,
      ymax = Inf
    )
  
  # Reference axis shown only in the upper-left facet.
  tick_values <- c(-plot_limit, 0, plot_limit)
  
  tick_labels <- c(paste0("-", formatC(
    plot_limit, format = "fg", digits = 2
  )),
  "0",
  formatC(plot_limit, format = "fg", digits = 2))
  
  axis_segment <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    x = 0,
    ymin = -plot_limit,
    ymax = plot_limit
  )
  
  tick_segments <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    y = tick_values
  )
  
  tick_text <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    y = tick_values,
    label = tick_labels
  )
  
  ggplot(plot_data, aes(x = trans, y = cc, fill = trans)) +
    geom_rect(
      data = background_negative,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "#3287a820"
    ) +
    geom_rect(
      data = background_positive,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "#a8324020"
    ) +
    geom_col(width = 0.9) +
    geom_segment(
      data = axis_segment,
      aes(
        x = x,
        xend = x,
        y = ymin,
        yend = ymax
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = tick_segments,
      aes(
        x = -0.2,
        xend = 0,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_text(
      data = tick_text,
      aes(x = -0.3, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 2.8
    ) +
    scale_fill_manual(
      values = transition_colors_s3,
      breaks = transition_levels_s3,
      labels = transition_labels_s3,
      drop = FALSE
    ) +
    scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0.015, 0.015))) +
    facet_grid(
      rows = vars(world1),
      cols = vars(world2),
      switch = "both",
      drop = FALSE
    ) +
    coord_cartesian(ylim = c(-plot_limit, plot_limit), clip = "off") +
    labs(x = "world 2", y = "world 1", fill = "hazard:") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      
      legend.position = "bottom",
      legend.box.spacing = unit(1, "mm"),
      
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.text.x.bottom = element_text(angle = 0),
      
      panel.spacing.x = unit(1.2, "mm"),
      panel.spacing.y = unit(1.2, "mm"),
      
      plot.margin = margin(5.5, 5.5, 5.5, 24)
    )
}

# -----------------------------------------------------------------------------
# Create and save both large figures
# -----------------------------------------------------------------------------

p_s3_returns <- make_s3_margin_plot("returns")

p_s3_noreturns <- make_s3_margin_plot("noreturns")

fig_s3_returns_file <- file.path("figs",
                                 paste0("dec_returns_margins_", interval_label, ".pdf"))

fig_s3_noreturns_file <- file.path("figs",
                                   paste0("dec_noreturns_margins_", interval_label, ".pdf"))

ggsave(
  filename = fig_s3_returns_file,
  plot = p_s3_returns,
  width = 30,
  height = 30,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  filename = fig_s3_noreturns_file,
  plot = p_s3_noreturns,
  width = 30,
  height = 30,
  units = "cm",
  device = cairo_pdf
)

p_s3_returns
p_s3_noreturns


# =============================================================================
# Supplementary Figure S4:
# Sullivan-comparable decomposition margins, all 16 worlds
# =============================================================================

# Combine the transition-specific contributions into the two components that
# would be reported by a Sullivan decomposition:
#   Mortality = h_hd + h_ud
#   Health    = h_hu + h_uh
#
# dec_s3_symmetric already contains both directions of every comparison, with
# the sign reversed correctly for the lower triangle.

component_levels_s4 <- c("Mortality", "Health")

component_colors_s4 <- c(
  Mortality = "#F8766D",
  Health    = "#00BFC4"
)

dec_s4_margins <-
  dec_s3_symmetric |>
  mutate(
    component = if_else(
      as.character(trans) %in% c("hd", "ud"),
      "Mortality",
      "Health"
    ),
    component = factor(component, levels = component_levels_s4)
  ) |>
  group_by(system, world1, world2, component) |>
  summarise(
    cc = sum(cc, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(world1 != world2)


# -----------------------------------------------------------------------------
# Function for one transition system
# -----------------------------------------------------------------------------

make_s4_sullivan_plot <- function(system_value) {
  plot_data <-
    dec_s4_margins |>
    filter(system == system_value)

  # Use a symmetric range, rounded upward to the next half-year, separately
  # for each system. This preserves visible variation within each large plot.
  plot_limit <- max(abs(plot_data$cc), na.rm = TRUE)
  plot_limit <- ceiling(plot_limit * 2) / 2

  if (!is.finite(plot_limit) || plot_limit == 0) {
    plot_limit <- 0.5
  }

  # All off-diagonal panels receive the positive/negative background. The
  # otherwise empty upper-left panel is retained for the reference y-axis.
  panel_grid <-
    expand_grid(
      world1 = factor(1:16, levels = 1:16),
      world2 = factor(1:16, levels = 1:16)
    ) |>
    filter(
      as.character(world1) != as.character(world2) |
        (as.character(world1) == "1" & as.character(world2) == "1")
    )

  background_negative <-
    panel_grid |>
    mutate(
      xmin = 0,
      xmax = 2.5,
      ymin = -Inf,
      ymax = 0
    )

  background_positive <-
    panel_grid |>
    mutate(
      xmin = 0,
      xmax = 2.5,
      ymin = 0,
      ymax = Inf
    )

  # Reference axis shown only in the upper-left facet.
  tick_values <- c(-plot_limit, 0, plot_limit)

  tick_labels <- c(
    paste0("-", formatC(plot_limit, format = "fg", digits = 2)),
    "0",
    formatC(plot_limit, format = "fg", digits = 2)
  )

  axis_segment <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    x = 0,
    ymin = -plot_limit,
    ymax = plot_limit
  )

  tick_segments <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    y = tick_values
  )

  tick_text <- tibble(
    world1 = factor(1, levels = 1:16),
    world2 = factor(1, levels = 1:16),
    y = tick_values,
    label = tick_labels
  )

  ggplot(plot_data, aes(x = component, y = cc, fill = component)) +
    geom_rect(
      data = background_negative,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "#3287a820"
    ) +
    geom_rect(
      data = background_positive,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "#a8324020"
    ) +
    geom_col(width = 0.78) +
    geom_segment(
      data = axis_segment,
      aes(
        x = x,
        xend = x,
        y = ymin,
        yend = ymax
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = tick_segments,
      aes(
        x = -0.2,
        xend = 0,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_text(
      data = tick_text,
      aes(x = -0.3, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 2.8
    ) +
    scale_fill_manual(
      values = component_colors_s4,
      breaks = component_levels_s4,
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = NULL,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    facet_grid(
      rows = vars(world1),
      cols = vars(world2),
      switch = "both",
      drop = FALSE
    ) +
    coord_cartesian(
      ylim = c(-plot_limit, plot_limit),
      clip = "off"
    ) +
    labs(
      x = "world 2",
      y = "world 1",
      fill = "component:"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),

      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),

      legend.position = "bottom",
      legend.box.spacing = unit(1, "mm"),

      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      strip.text.x.bottom = element_text(angle = 0),

      panel.spacing.x = unit(1.2, "mm"),
      panel.spacing.y = unit(1.2, "mm"),

      plot.margin = margin(5.5, 5.5, 5.5, 24)
    )
}


# -----------------------------------------------------------------------------
# Create and save both large figures
# -----------------------------------------------------------------------------

p_s4_sullivan_returns <- make_s4_sullivan_plot("returns")

p_s4_sullivan_noreturns <- make_s4_sullivan_plot("noreturns")

fig_s4_returns_file <-
  file.path(
    "figs",
    paste0("dec_sullivan_returns_margins_", interval_label, ".pdf")
  )

fig_s4_noreturns_file <-
  file.path(
    "figs",
    paste0("dec_sullivan_noreturns_margins_", interval_label, ".pdf")
  )

ggsave(
  filename = fig_s4_returns_file,
  plot = p_s4_sullivan_returns,
  width = 30,
  height = 30,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  filename = fig_s4_noreturns_file,
  plot = p_s4_sullivan_noreturns,
  width = 30,
  height = 30,
  units = "cm",
  device = cairo_pdf
)

p_s4_sullivan_returns
p_s4_sullivan_noreturns
