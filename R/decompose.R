
suppressPackageStartupMessages(library(tidyverse))

haz <- read_csv("data/worlds_hazards_annual.csv.gz")
haz

worlds <- haz$world |> unique() |> sort()
world_pairings <- expand_grid(world1 = worlds, world2 = worlds) |> 
  filter(world1 < world2)

world_pairings <- bind_rows(world_pairings |> mutate(system = "returns"),
                            world_pairings |> mutate(system = "noreturns"))

haz2 <- haz|> 
  rename(hazard2 = hazard,
         world2 = world)
haz <- haz |> 
  rename(hazard1 = hazard,
         world1 = world)

dec_data <- 
  full_join(haz, world_pairings, by = join_by(world1, system),relationship = "many-to-many") |> 
  left_join(haz2, by = join_by(age,trans,world2, system)) |> 
  select(system,world1, world2,trans,age,hazard1,hazard2) |> 
  arrange(system,world1, world2,trans,age)

