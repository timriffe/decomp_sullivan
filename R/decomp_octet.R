# this script currently sets up how a single decomposition might work.
# it would


library(tidyverse)
library(DemoDecomp)
source("R/simulation_functions.R")

octet <- read_csv("data/haz_octet.csv.gz")
head(octet)

mslt <-
octet |> 
  filter(world == 5) |> 
  haz_to_probs(age = 0:100, age_int = 1) |> 
  calculate_lt( init =  c(H=1,U=0), age_int = 1)
  
mslt |> 
  summarize(hle = sum((1-prevalence)*lx),
            hle_ms = sum(lh),
            d = hle - hle_ms)
head(mslt)

noreturns <- octet |> filter(system == "noreturn")


haz_to_hle <- function(haz_vec, trans, age, age_int = 1, init = c(H=1,U=0)){
  data <- tibble(trans = trans, age = age, hazard = haz_vec)
  lh <- 
    data |> 
    haz_to_probs(age = age, age_int = age_int) |> 
    calculate_lt( init = init) |> 
    pull(lh)
  sum((lh + c(lh[-1],0))/2)
}

returns |> 
  group_by(world) |> 
  summarize(hle = haz_to_hle(haz_vec = hazard,
                             age = age,
                             trans = trans,
                             age_int = 1,
                             init = c(H=1,U=0))) |> 
  pull(hle)



worlds <- unique(returns$world)
worlds

dec_sets <- combn(worlds,2)
as.symbol()
for (i in 1:ncol(dec_sets)){
  wrldsi <- dec_sets[, i] |> as.character()
  
  dec_i <- returns |> 
    filter(world %in% wrldsi) |> 
    pivot_wider(names_from = world, values_from = hazard) |> 
    mutate(cc = horiuchi(func = haz_to_hle,
                         pars1 = get(wrldsi[1]),
                         pars2 = get(wrldsi[2]),
                         trans = trans,
                         age = age,
                         age_int = 1,
                         init = c(H=1,U=0),
                         N = 100))
}

dec_i |> 
  ggplot(aes(x =age, y = cc, color = trans)) +
  geom_line()

dec_i$cc |> sum()


dec_i |> 
  group_by(sign(cc)) |> 
  summarize(cc = sum(cc))

