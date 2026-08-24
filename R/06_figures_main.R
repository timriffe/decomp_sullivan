# fig code might need reordering

library(tidyverse)
library(patchwork)

# ============================================================
# Main-manuscript figure settings
# ============================================================

# Use 1 for annual decomposition or 1/12 for monthly
# decomposition summed back to completed single ages.
age_int <- 1 / 12

# Four representative worlds for the manuscript figures.
main_worlds <- c(2, 6, 11, 14)

# Shared transition order and display labels.
# Shared transition order and plotmath labels.
transition_levels <- c("hu", "uh", "hd", "ud")

transition_labels <- c(
  hu = "h[hu]",
  uh = "h[uh]",
  hd = "h[hd]",
  ud = "h[ud]"
)

transition_expressions <- parse(
  text = unname(transition_labels)
)

# Shared output dimensions for both manuscript figures.
fig_width  <- 20
fig_height <- 15

# Relative width of the blank space between Returns and No returns.
system_gap <- 0.08

dir.create("figs", showWarnings = FALSE, recursive = TRUE)

age_label <- if_else(age_int == 1, "annual", "monthly")

fig1_name <- paste0(
  "figs/dec_margins_",
  age_label,
  "_main.pdf"
)

fig2_name <- paste0(
  "figs/dec_sullivan_margins_",
  age_label,
  "_main.pdf"
)

fig3_name <- paste0(
  "figs/world_hazards_",
  age_label,
  "_main.pdf"
)


# ============================================================
# Read and prepare decomposition results
# ============================================================

if (age_int == 1) {
  
  dec <- read_csv(
    "data/worlds_decomp_annual.csv.gz",
    show_col_types = FALSE
  )
  haz <- read_csv(
    "data/worlds_hazards_annual.csv.gz",
    show_col_types = FALSE
  )
}

if (age_int == 1 / 12) {
  
  dec <- read_csv(
    "data/worlds_decomp_monthly.csv.gz",
    show_col_types = FALSE
  ) |>
    # Sum monthly contributions back to completed single ages.
    mutate(age = age - age %% 1) |>
    summarize(
      cc = sum(cc, na.rm = TRUE),
      .by = c(system, world1, world2, trans, age)
    )
  haz <- read_csv(
    "data/worlds_hazards_monthly.csv.gz",
    show_col_types = FALSE
  )
}

dec <- dec |>
  filter(
    world1 %in% main_worlds,
    world2 %in% main_worlds,
    system %in% c("returns", "noreturns")
  ) |>
  mutate(
    system = factor(
      system,
      levels = c("returns", "noreturns"),
      labels = c("Returns", "No returns")
    ),
    world1 = factor(world1, levels = main_worlds),
    world2 = factor(world2, levels = main_worlds),
    trans = factor(trans, levels = transition_levels)
  )


# Add reverse comparisons for the lower triangle.
dec_symmetric <- dec |>
  mutate(
    world_hold = world1,
    world1 = world2,
    world2 = world_hold,
    cc = -cc
  ) |>
  select(-world_hold) |>
  bind_rows(dec)


# ============================================================
# Figure 1: transition-specific margins
# ============================================================

dec_fig_margins <- dec_symmetric |>
  summarize(
    cc = sum(cc, na.rm = TRUE),
    .by = c(system, world1, world2, trans)
  ) |>
  # Diagonal panels carry no comparison.
  filter(world1 != world2)


# Common symmetric y range for Returns and No returns.
lim_transition <- max(abs(dec_fig_margins$cc), na.rm = TRUE)

if (!is.finite(lim_transition) || lim_transition == 0) {
  lim_transition <- 0.5
} else {
  lim_transition <- ceiling(lim_transition * 1.05 * 2) / 2
}


# Background rectangles for the non-diagonal panels.
bg_transition <- expand_grid(
  world1 = factor(main_worlds, levels = main_worlds),
  world2 = factor(main_worlds, levels = main_worlds)
) |>
  filter(world1 != world2)

bg_transition_neg <- bg_transition |>
  mutate(
    xmin = 0,
    xmax = 4.5,
    ymin = -lim_transition,
    ymax = 0
  )

bg_transition_pos <- bg_transition |>
  mutate(
    xmin = 0,
    xmax = 4.5,
    ymin = 0,
    ymax = lim_transition
  )


# ----------------------------
# Transition margins: Returns
# ----------------------------

p_transition_returns <-
  dec_fig_margins |>
  filter(system == "Returns") |>
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  
  geom_rect(
    data = bg_transition_neg,
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
    data = bg_transition_pos,
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
  
  scale_fill_discrete(
    breaks = transition_levels,
    labels = transition_expressions,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(-lim_transition, lim_transition),
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  facet_grid(
    rows = vars(world1),
    cols = vars(world2),
    switch = "both",
    drop = FALSE
  ) +
  
  labs(
    title = "Returns",
    x = "world 2",
    y = "world 1",
    fill = "hazard:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    panel.spacing.x = unit(2.2, "mm"),
    panel.spacing.y = unit(2.2, "mm"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.text.x.bottom = element_text(angle = 0),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )


# -------------------------------
# Transition margins: No returns
# -------------------------------

p_transition_noreturns <-
  dec_fig_margins |>
  filter(system == "No returns") |>
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  
  geom_rect(
    data = bg_transition_neg,
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
    data = bg_transition_pos,
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
  scale_fill_discrete(
    breaks = transition_levels,
    labels = transition_expressions,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(-lim_transition, lim_transition),
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  facet_grid(
    rows = vars(world1),
    cols = vars(world2),
    switch = "both",
    drop = FALSE
  ) +
  
  labs(
    title = "No returns",
    x = "world 2",
    y = NULL,
    fill = "hazard:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # The left block already supplies the shared y-axis labels.
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    panel.spacing.x = unit(2.2, "mm"),
    panel.spacing.y = unit(2.2, "mm"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    strip.text.x.bottom = element_text(angle = 0),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 2)
  )


p_transition <-
  p_transition_returns +
  plot_spacer() +
  p_transition_noreturns +
  plot_layout(
    widths = c(1, system_gap, 1),
    guides = "collect"
  ) &
  theme(legend.position = "bottom")


p_transition

ggsave(
  filename = fig1_name,
  plot = p_transition,
  width = fig_width,
  height = fig_height,
  units = "cm",
  device = cairo_pdf
)


# ============================================================
# Figure 2: Sullivan-comparable mortality and health margins
# ============================================================

dec_fig_margins_sullivan <- dec_symmetric |>
  mutate(
    component = if_else(
      trans %in% c("hd", "ud"),
      "mortality",
      "health"
    ),
    component = factor(
      component,
      levels = c("mortality", "health"),
      labels = c("Mortality", "Health")
    )
  ) |>
  summarize(
    cc = sum(cc, na.rm = TRUE),
    .by = c(system, world1, world2, component)
  ) |>
  filter(world1 != world2)


# Common symmetric y range for Returns and No returns.
lim_sullivan <- max(abs(dec_fig_margins_sullivan$cc), na.rm = TRUE)

if (!is.finite(lim_sullivan) || lim_sullivan == 0) {
  lim_sullivan <- 0.5
} else {
  lim_sullivan <- ceiling(lim_sullivan * 1.05 * 2) / 2
}


# Background rectangles for the non-diagonal panels.
bg_sullivan <- expand_grid(
  world1 = factor(main_worlds, levels = main_worlds),
  world2 = factor(main_worlds, levels = main_worlds)
) |>
  filter(world1 != world2)

bg_sullivan_neg <- bg_sullivan |>
  mutate(
    xmin = 0,
    xmax = 2.5,
    ymin = -lim_sullivan,
    ymax = 0
  )

bg_sullivan_pos <- bg_sullivan |>
  mutate(
    xmin = 0,
    xmax = 2.5,
    ymin = 0,
    ymax = lim_sullivan
  )


# --------------------------
# Sullivan margins: Returns
# --------------------------

p_sullivan_returns <-
  dec_fig_margins_sullivan |>
  filter(system == "Returns") |>
  ggplot(aes(x = component, y = cc, fill = component)) +
  
  geom_rect(
    data = bg_sullivan_neg,
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
    data = bg_sullivan_pos,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "#a8324020"
  ) +
  
  geom_col(width = 0.68) +
  
  scale_y_continuous(
    limits = c(-lim_sullivan, lim_sullivan),
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  facet_grid(
    rows = vars(world1),
    cols = vars(world2),
    switch = "both",
    drop = FALSE
  ) +
  
  labs(
    title = "Returns",
    x = "world 2",
    y = "world 1",
    fill = "component:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    panel.spacing.x = unit(2.2, "mm"),
    panel.spacing.y = unit(2.2, "mm"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.text.x.bottom = element_text(angle = 0),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )


# -----------------------------
# Sullivan margins: No returns
# -----------------------------

p_sullivan_noreturns <-
  dec_fig_margins_sullivan |>
  filter(system == "No returns") |>
  ggplot(aes(x = component, y = cc, fill = component)) +
  
  geom_rect(
    data = bg_sullivan_neg,
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
    data = bg_sullivan_pos,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "#a8324020"
  ) +
  
  geom_col(width = 0.68) +
  
  scale_y_continuous(
    limits = c(-lim_sullivan, lim_sullivan),
    breaks = scales::breaks_pretty(n = 3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  facet_grid(
    rows = vars(world1),
    cols = vars(world2),
    switch = "both",
    drop = FALSE
  ) +
  
  labs(
    title = "No returns",
    x = "world 2",
    y = NULL,
    fill = "component:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    panel.spacing.x = unit(2.2, "mm"),
    panel.spacing.y = unit(2.2, "mm"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    strip.text.x.bottom = element_text(angle = 0),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 2)
  )


p_sullivan <-
  p_sullivan_returns +
  plot_spacer() +
  p_sullivan_noreturns +
  plot_layout(
    widths = c(1, system_gap, 1),
    guides = "collect"
  ) &
  theme(legend.position = "bottom")


p_sullivan

ggsave(
  filename = fig2_name,
  plot = p_sullivan,
  width = fig_width,
  height = fig_height,
  units = "cm",
  device = cairo_pdf
)


# ============================================================
# Figure 3: age-specific hazards in the four selected worlds
# ============================================================

# Worlds carrying the same number do not represent the same scenario
# in the Returns and No returns systems. We therefore make a separate
# plot block and legend for each system rather than using one shared
# world-colour scale across systems.

haz <- haz |>
  filter(
    world %in% main_worlds,
    system %in% c("returns", "noreturns")
  ) |>
  mutate(
    trans = factor(
      trans,
      levels = transition_levels,
      labels = unname(transition_labels)
    ),
    world = factor(
      world,
      levels = main_worlds
    )
  )

# ----------------------------
# Hazard profiles: Returns
# ----------------------------

p_haz_returns <-
  haz |>
  filter(system == "returns") |>
  ggplot(aes(x = age, y = hazard, colour = world)) +
  
  geom_line(linewidth = 0.65) +
  
  facet_grid(
    rows = vars(trans),
    scales = "free_y",
    drop = FALSE,
    switch = "y",
    labeller = labeller(
      trans = label_parsed
    )
  ) +
  
  scale_colour_brewer(
    palette = "Dark2",
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(50, 100, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_log10(
    labels = function(x) {
      sub(
        "\\.?0+$",
        "",
        format(x, scientific = FALSE, trim = TRUE)
      )
    }
  ) +
  
  labs(
    title = "Returns",
    x = "Age",
    y = "Hazard",
    colour = "Returns world:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    axis.text.y = element_text(hjust = 0),
    
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(2.2, "mm"),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )


# -------------------------------
# Hazard profiles: No returns
# -------------------------------

p_haz_noreturns <-
  haz |>
  filter(system == "noreturns") |>
  ggplot(aes(x = age, y = hazard, colour = world)) +
  
  geom_line(linewidth = 0.65) +
  
  facet_grid(
    rows = vars(trans),
    scales = "free_y",
    drop = FALSE,
    switch = "y"
  ) +
  
  # A different palette and legend title avoid implying that, for example,
  # Returns world 1 and No returns world 1 are the same scenario.
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = seq(50, 100, by = 10),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_y_log10(
    labels = function(x) {
      sub(
        "\\.?0+$",
        "",
        format(x, scientific = FALSE, trim = TRUE)
      )
    }
  ) +
  
  labs(
    title = "No returns",
    x = "Age",
    y = "",
    colour = "No returns world:"
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(2.2, "mm"),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm"),
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(5.5, 5.5, 5.5, 2)
  )


p_hazards <-
  p_haz_returns +
  plot_spacer() +
  p_haz_noreturns +
  plot_layout(
    widths = c(1, system_gap, 1),
    guides = "keep"
  )

p_hazards

ggsave(
  filename = fig3_name,
  plot = p_hazards,
  width = fig_width,
  height = fig_height,
  units = "cm",
  device = cairo_pdf
)