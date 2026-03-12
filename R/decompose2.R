suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(furrr)
  library(future)
  library(DemoDecomp)
  library(progressr)
})

source("R/decompose_functions2.R")

test <- FALSE

# globals
N <- 20
age_int <- 1

if (age_int == 1) {
  haz <- read_csv("data/worlds_hazards_annual.csv.gz", show_col_types = FALSE)
}
if (age_int == (1 / 12)) {
  haz <- read_csv("data/worlds_hazards_monthly.csv.gz", show_col_types = FALSE)
}

if (test) {
  haz <- haz |> filter(world < 4)
  N <- 2
}

init <- c(h = 1, u = 0)
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
haz2 <- haz |> rename(hazard2 = hazard, world2 = world)
haz1 <- haz |> rename(hazard1 = hazard, world1 = world)

dec_data <- full_join(
  haz1,
  world_pairings,
  by = join_by(world1, system),
  relationship = "many-to-many"
) |>
  left_join(haz2, by = join_by(age, trans, world2, system)) |>
  select(system, world1, world2, trans, age, hazard1, hazard2) |>
  arrange(system, world1, world2, trans, age) |>
  filter(world1 < world2)

rm(list = c("haz", "haz1", "haz2"))
gc()

# ready
#
# Decomp in parallel
plan(multisession, workers = n_cores)

# Split into list of tibbles, one per group
dec_pairs <- dec_data |>
  group_by(world1, world2, system) |>
  group_split(.keep = TRUE)

rm(dec_data)
gc()

# Helper that runs ONE group (one custom horiuchi call) and returns that group's rows + cc column
do_decomp <- function(g, N, age_int, init, p = NULL) {
  sys <- dplyr::first(g$system)
  
  out <- g |>
    dplyr::mutate(
      cc = horiuchi_haz_cached(
        pars1 = hazard1,
        pars2 = hazard2,
        N = N,
        trans = trans,
        age = age,
        age_int = age_int,
        init = init,
        system = sys
      )
    )
  
  # for progress bar
  if (!is.null(p)) {
    p()
  }
  
  return(out)
}

handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = dec_pairs)
  
  # Parallel map over groups; then bind back into one tibble
  dec_out <- future_map(
    dec_pairs,
    do_decomp,
    N = N,
    age_int = age_int,
    init = init,
    p = p,
    .options = furrr::furrr_options(seed = TRUE)
  ) |>
    bind_rows()
})

plan(sequential)

# clean up
# save out
if (age_int == 1) {
  write_csv(dec_out, "data/worlds_decomp_annual.csv.gz")
}
if (age_int == 1 / 12) {
  write_csv(dec_out, "data/worlds_decomp_monthly.csv.gz")
}

# end
dec_out |>
  group_by(world1, world2, system) |>
  summarize(cc = sum(cc))