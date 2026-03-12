

library(tidyverse)

age_int <- 1/12
system <- "returns"

fig1_name <- paste0("figs/dec_",system,"_margins_",ifelse(age_int==1,"annual.pdf","monthly.pdf"))
fig2_name <- paste0("figs/dec_",system,"_age_",ifelse(age_int==1,"annual.pdf","monthly.pdf"))

if (age_int == 1){
  dec <- read_csv("data/worlds_decomp_annual.csv.gz")
}
if (age_int == 1/12){
  dec <- read_csv("data/worlds_decomp_monthly.csv.gz")
  dec <- dec |> 
    mutate(age = age - age %% 1) |> 
    group_by(system, world1, world2, trans, age) |> 
    summarize(cc = sum(cc), .groups = "drop")
}

dec |> 
  filter(system == "returns") |> 
  ggplot(aes(x = age, y = cc, color = trans)) +
  geom_line() +
  facet_grid(vars(world1),vars(world2))

dec |> 
  filter(system == "noreturns") |> 
  ggplot(aes(x = age, y = cc, color = trans)) +
  geom_line() +
  facet_grid(vars(world1),vars(world2))



dec_fig_margins <-
  dec |>
  mutate(world_hold = world1,
         world1 = world2,
         world2 = world_hold,
         cc = -cc) |>
  select(-world_hold) |>
  bind_rows(dec)  |> 
  filter(system == "returns") |>
  summarize(cc = sum(cc),
            .by = c(world1, world2, trans))

tick_dat <- data.frame(
  world1 = 1,
  world2 = 1,
  y = c(-3, 0, 3),
  lab = c("-3", "0", "3")
)

tick_seg <- data.frame(
  world1 = 1,
  world2 = 1,
  y = c(-3, 0, 3)
)

axis_seg <- data.frame(
  world1 = 1,
  world2 = 1,
  x = 0
)
bg_dat <- expand_grid(
  world1 = sort(unique(dec$world1)),
  world2 = sort(unique(dec$world2))
) |>
  filter(world1 != world2 | (world1 == 1 & world2 == 1))

bg_bot <- bg_dat |>
  mutate(ymin = -Inf, ymax = 0,
         xmin = 0, xmax = 4.5,
         fill_col = "#3287a820")

bg_top <- bg_dat |>
  mutate(ymin = 0, ymax = Inf,
         xmin = 0, xmax = 4.5,
         fill_col = "#a8324020")
bg_panels <- dec_fig_margins |>
  distinct(world1, world2) |>
  full_join(
    expand_grid(
      world1 = sort(unique(dec_fig_margins$world1)),
      world2 = sort(unique(dec_fig_margins$world2))
    ),
    by = c("world1", "world2")
  ) |>
  distinct(world1, world2) |>
  filter(world1 != world2 | (world1 == 1 & world2 == 1))

bg_neg <- bg_panels |>
  mutate(xmin = 0, xmax = 4.5, ymin = -Inf, ymax = 0)

bg_pos <- bg_panels |>
  mutate(xmin = 0, xmax = 4.5, ymin = 0, ymax = Inf)

p<-
dec_fig_margins |>
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  
  geom_rect(
    data = bg_neg,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#3287a820"
  ) +
  geom_rect(
    data = bg_pos,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#a8324020"
  ) +
  
  geom_col() +
  
  geom_segment(
    data = axis_seg,
    aes(x = x, xend = x, y = -3.1, yend = 3.1),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = tick_seg,
    aes(x = -.2, xend = 0, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tick_dat,
    aes(x = -.3, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.8
  ) +
  
  scale_y_continuous(breaks = NULL) +
  facet_grid(vars(world1), vars(world2), switch = "both") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 5.5, 22),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x.bottom = element_text(angle = 0)
  ) +
  labs(x = "world 2", y = "world 1", fill = "hazard:")


ggsave(p,file=fig1_name,width=30,height=30,units="cm")



# --------------------------------------
# age profile panel:
dec_all <-
  dec |>
  mutate(world_hold = world1,
         world1 = world2,
         world2 = world_hold,
         cc = -cc) |>
  select(-world_hold) |>
  bind_rows(dec) 


bg_panels <- expand_grid(
  world1 = sort(unique(dec_all$world1)),
  world2 = sort(unique(dec_all$world2))
) |>
  filter(world1 != world2 | (world1 == 1 & world2 == 1))

bg_neg <- bg_panels |>
  mutate(xmin = 50, xmax = 100, ymin = -0.1, ymax = 0)

bg_pos <- bg_panels |>
  mutate(xmin = 50, xmax = 100, ymin = 0, ymax = 0.1)

# fake y-axis for upper-left facet
tick_dat_y <- data.frame(
  world1 = 1,
  world2 = 1,
  y = c(-0.1, 0, 0.1),
  lab = c("-0.1", "0", "0.1")
)

tick_seg_y <- data.frame(
  world1 = 1,
  world2 = 1,
  y = c(-0.1, 0, 0.1)
)

axis_seg_y <- data.frame(
  world1 = 1,
  world2 = 1,
  x = 50
)

# fake x-axis for upper-left facet
tick_dat_x <- data.frame(
  world1 = 1,
  world2 = 1,
  x = c(50, 75, 100),
  lab = c("50", "Age", "100")
)

tick_seg_x <- data.frame(
  world1 = 1,
  world2 = 1,
  x = c(50, 100),
  y0 = c(-0.1, -0.1),
  y1 = c(-0.106, -0.05)
)

axis_seg_x <- data.frame(
  world1 = 1,
  world2 = 1,
  y = -0.1
)

p<-
dec_all |>
  filter(system=="returns") |> 
  ggplot(aes(x = age, y = cc, color = trans)) +
  
  geom_rect(
    data = bg_neg,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#3287a820"
  ) +
  geom_rect(
    data = bg_pos,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#a8324020"
  ) +
  
  geom_line() +
  
  # fake y-axis only in upper-left facet
  geom_segment(
    data = axis_seg_y,
    aes(x = x, xend = x, y = -0.1, yend = 0.1),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = tick_seg_y,
    aes(x = 49.2, xend = 50, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tick_dat_y,
    aes(x = 48.9, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.8
  ) +
  
  # fake x-axis only in upper-left facet
  geom_segment(
    data = tick_seg_x,
    aes(x = x, xend = x, y = y0, yend = y1),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = tick_dat_x,
    aes(x = x, y = -0.112, label = lab),
    inherit.aes = FALSE,
    vjust = 1,
    size = 2.8
  ) +
  
  scale_x_continuous(
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  facet_grid(vars(world1), vars(world2), switch = "both") +
  coord_cartesian(
    xlim = c(48, 100),
    ylim = c(-0.12, 0.1),
    clip = "off"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(12, 5.5, 18, 30),
    strip.placement = "outside",
    strip.background = element_blank(),
    
    strip.text.y.left = element_text(angle = 0),
    strip.text.x.bottom = element_text(angle = 0)
  ) +
  labs(x = "world 2", y = "world 1", color = "hazard:")


ggsave(p,file=fig2_name,width=30,height=30,units="cm")
