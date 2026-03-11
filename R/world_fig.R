

library(tidyverse)
dec <- read_csv("data/worlds_decomp_annual.csv.gz")

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



dec |> 
  filter(system == "returns") |> 
  group_by(world1,world2,trans) |> 
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  geom_col() +
  facet_grid(vars(world1),vars(world2))+
  theme_minimal()

# install.packages("ggh4x")
library(ggh4x)
dec |> 
  mutate(world_hold = world1,
         world1 = world2,
         world2 = world_hold,
         cc = -cc) |> 
  select(-world_hold) |> 
  bind_rows(dec) |> 
  filter(system == "returns") |> 
  group_by(world1,world2,trans) |> 
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  geom_rect(data = data.frame(x = 0, y = 0,trans=NA,cc=0),aes(xmin=0,xmax=4.5,ymin=-Inf,ymax=0),
            fill="#3287a820")+
  geom_rect(data = data.frame(x = 0, y = 0,trans=NA,cc=0),aes(xmin=0,xmax=4.5,ymin=0,ymax=Inf),
            fill="#a8324020")+
  geom_col() +
  scale_y_continuous(breaks = c(-3,0,3))+
  facet_grid(vars(world1),vars(world2))+
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position="bottom"
        ) +
  labs(x="world 2", y = "world 1",fill="hazard:")
library(dplyr)
library(ggplot2)

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
dec |>
  mutate(world_hold = world1,
         world1 = world2,
         world2 = world_hold,
         cc = -cc) |>
  select(-world_hold) |>
  bind_rows(dec) |>
  filter(system == "returns") |>
  group_by(world1, world2, trans) |>
  ggplot(aes(x = trans, y = cc, fill = trans)) +
  
  annotate("rect", xmin = 0, xmax = 4.5, ymin = -Inf, ymax = 0,
           fill = "#3287a820") +
  annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = Inf,
           fill = "#a8324020") +
  
  geom_col() +
  
  # fake y-axis only in upper-left facet
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
  facet_grid(vars(world1), vars(world2),switch = "both") +
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
    strip.background = element_blank()
  ) +
  labs(x = "world 2", y = "world 1", fill = "hazard:")
