
# Not sure, but methinks the sd of healthy years, calculated per ms_dist functions, will differ between worlds as well, and this might be worth pointing out.


suppressPackageStartupMessages({
  library(tidyverse)
  library(collapse)
  library(tidyfast)
})
source("https://raw.githubusercontent.com/timriffe/ms_dist/refs/heads/master/code/01_functions.R")

# calc_dxh() wants to see a p_tibble. It's also premised on single ages. Let's generalize a bit?
calc_dxh <- function(p_tibble, age_int = 1, init = NULL){
  ages <- p_tibble$age |> unique() |> sort()
  aa   <- c(ages,max(ages)+age_int)
  min_age <- min(aa)
  hh   <- aa - min(aa)
  
  if (is.null(init)){
  init <- p_tibble |> 
    filter(age == min_age) |> 
    init_constant()
  }
  d3 <- expand.grid(h = hh,
                    age = aa,
                    current_state = c("H","U"),
                    l = 0,
                    stringsAsFactors = FALSE) |>
    fsubset(h <= (age - min_age)) |>
    fmutate(l = data.table::fcase(
      age == min_age & current_state == "H" , init["H"],
      age == min_age & current_state == "U" , init["U"],
      default = 0))
  for (a in ages) {
    
    d3n         <- fsubset(d3, age == a)
    possible_ds <- unique(d3n$h) |> sort()
    dt          <- fsubset(p_tibble,age == a)
    
    HH <- dt$HH
    HU <- dt$HU
    UU <- dt$UU
    UH <- dt$UH
    
    for (d in possible_ds) {
      
      d3nd <- fsubset(d3n, h == d)
      lxhH <- fsubset(d3nd, current_state == "H")$l
      lxhU <- fsubset(d3nd, current_state == "U")$l
      
      d3 <- d3 |>
        fmutate(l = fcase(
          abs(age - (a + age_int)) < 1e-5 & h == d + age_int & current_state == "H", l + lxhH * HH + lxhU * UH,
          abs(age - (a + age_int)) < 1e-5 & h == d & current_state == "U", l + lxhU * UU + lxhH * HU,
          rep_len(TRUE, length(l)), l))
    }
  }
  d_out <- p_tibble |>
    select(age, HD, UD) |>
    rename(H = HD, 
           U = UD) |>
    pivot_longer(c(H, U), 
                 names_to  = "current_state", 
                 values_to = "qx") |> 
    right_join(d3, by = c("age", "current_state"), multiple = "all") |>
    mutate(qx  = ifelse(age == max(age), 1, qx),
           dxs = qx * l,
           x   = age - min_age,
           u   = x - h) |> 
    select(current_state, age, x, h, u, lxsc = l, dxsc = dxs)
  d_out
}
age_int <- 1
if (age_int == 1){
  worlds <- read_csv("data/worlds_mslt_annual.csv.gz")  
} 
if (age_int == 1/12){
  worlds <- read_csv("data/worlds_mslt_monthly.csv.gz")
}

# if age_int is 1/12, this can take a long time.
dxh <-
  worlds |> 
  select(system, world, age, HH, HU, HD, UH, UU, UD) |> 
  group_by(system, world) |> 
  do(calc_dxh(p_tibble = .data, age_int = age_int, init = c(H=1,U=0))) |> 
  ungroup()
dxh
# Problem: something must be different about how lxs is calculated. yeesh. Because we need all expectancies to be exactly equal. So something is off.
variances <-
  dxh |>
  summarize(
    le = sum((x + .5) * dxsc),
    hle = sum((h + .5) * dxsc),
    ule = sum((u + .5) * dxsc),
    vle    = sum(((x+.5) - le) ^ 2 * dxsc),
    vh     = sum(((h+.5) - hle) ^ 2 * dxsc),
    vu     = sum(((u+.5) - ule) ^ 2 * dxsc),
    cov_hu = sum(((h+.5) - hle) * ((u+.5) - ule) * dxsc),
    .by = c(system, world)) |>
  mutate(vle_check = vh + vu + 2 * cov_hu) 
variances |> View()


# Covariance between years lived in H and years lived in U
# is inconsistent between worlds.
variances |>
  ggplot(aes(x = world, y = cov_hu)) +
  geom_point() +
  facet_wrap(~system)

# We see:
# 1) total sd of lifespan is identical between worlds.
# 2) sd of total healthy years lived is identical in no-returns worlds only. It differs over the returns worlds
# 3) sd of unhealthy years is inconsistent between all worlds.
variances |> 
  mutate(sdx = sqrt(vle),
         sdh = sqrt(vh),
         sdu = sqrt(vu)) |> 
  select(system, world, sdx, sdh, sdu) |> 
  pivot_longer(sdx:sdu, names_to = "state", values_to = "sd") |> 
  ggplot(aes(x = world, y = sd, color=state)) +
  geom_point() +
  facet_wrap(~system)

