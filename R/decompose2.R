source("R/decompose_functions2.R")

test=FALSE
# globals
N     <- 20
age_int = 1
if (age_int == 1) {
  haz <- read_csv("data/worlds_hazards_annual.csv.gz", show_col_types = FALSE)
}
if (age_int == (1 / 12)) {
  haz <- read_csv("data/worlds_hazards_monthly.csv.gz", show_col_types = FALSE)
}

init  <- c(h = 1, u = 0)
n_cores <- parallel::detectCores()

# create world combinations
worlds <- haz$world |> unique() |> sort()
world_pairings <- expand_grid(world1 = worlds, world2 = worlds) |>
  filter(world1 < world2)

world_pairings <- bind_rows(
  world_pairings |> mutate(system = "returns"),
  world_pairings |> mutate(system = "noreturns")
)

# reformat data to have world pairings side by side
haz2 <- haz |>
  rename(hazard2 = hazard, world2 = world)
haz1 <- haz |>
  rename(hazard1 = hazard, world1 = world)

dec_data <-
  full_join(haz1,
            world_pairings,
            by = join_by(world1, system),
            relationship = "many-to-many") |>
  left_join(haz2, by = join_by(age, trans, world2, system)) |>
  select(system, world1, world2, trans, age, hazard1, hazard2) |>
  arrange(system, world1, world2, trans, age) |>
  filter(world1 < world2)

rm(list=c("haz", "haz1", "haz2"))
gc()
#
plan(multisession, workers = n_cores)

# Split into list of tibbles, one per group
dec_pairs <- dec_data |>
  group_by(world1, world2, system) |>
  group_split(.keep = TRUE)

rm(dec_data)
gc()