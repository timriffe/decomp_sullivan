
# install.packages("ipumsr")
library(ipumsr)
library(tidyverse)
library(mgcv)

# R.utils::gunzip("nhis_00002.dat.gz", remove = FALSE)
ddi         <- read_ipums_ddi("nhis_00002.xml")
data        <- read_ipums_micro(ddi)
names(data) <- tolower(names(data))
# ipums_view(ddi)

# weight
# PERWEIGHT
# mortality technically has its own corrected weight
# but it is present only till year 2015 combining 2 weights would be strange,
# so I resolve to using only original weight

# some data cleaning
dt <- data %>% 
  dplyr::select(year, id = nhispid, wt = perweight, age:hyp2time, 
                mortstat, mortdodq, mortdody, wt_m = mortwt, mortwtsa) %>% 
  as_factor()

# wt is roughly equal to us population size. so we can calculate rates directly
# weight for 2019 is missing so we omit this year
dt %>% 
  group_by(year) %>% 
  summarise(pop = sum(wt, na.rm = TRUE))

# Top code for 85 years or older (1963-1968 and 1997-forward)
# The open age group is 85+ unfortunately. 
table(dt$age)

# unique id correspondence. data is cross-sectional. no follow up possible
length(unique(dt$id))
nrow(dt)

# we can omit these person no important information here
dt %>% 
  filter(age %in% c("Unknown-don't know", "Unknown-refused"))

# proceed with cleaning
# wt_m adjusted to exclude NIU population.
# I think this is what we should use for analysis, but it only covers years < 2015
# data is missing for 2019
dt <- dt %>%
  dplyr::select(-c(mortdodq, mortwtsa)) %>% 
  filter(!age %in% c("Unknown-don't know", "Unknown-refused")) %>% 
  mutate(age = ifelse(str_detect(age, "85"), "85", age),
         age = as.numeric(age)) %>%
  # this period has almost no information on dementia
  dplyr::select(-dementiaev) %>% 
  # remove year with missing weights
  filter(year < 2019) %>% 
  mutate(year = as.numeric(as.character(year)))

# overall population
pop <- dt %>%
  #filter(mortstat != "NIU") |> 
  group_by(year, sex, age) %>%
  summarise(pop = sum(wt), .groups = "drop") %>% 
  pivot_wider(names_from = sex,
              values_from = pop) %>% 
  mutate(Total = Male + Female) %>% 
  pivot_longer(c(Male:Total),
               names_to = "sex",
               values_to = "pop")

# now lets calculate deaths
deaths <- dt %>%
  filter(mortstat == "Assumed deceased") %>%
  dplyr::select(wt, id, age, sex:epilepsyev, year = mortdody) %>%
  pivot_longer(-c(id:sex, year, wt),
               names_to = "cause",
               values_to = "status") %>%
  filter(status %in% c("No", "Yes")) %>% 
  group_by(year, sex, cause, status, age) %>% 
  summarise(death = sum(wt), .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year))) %>%
  filter(year < 2019) %>%
  pivot_wider(names_from = sex,
              values_from = death) %>%
  replace(is.na(.), 0) %>% 
  mutate(Total = Male + Female) %>% 
  pivot_longer(c(Male:Total),
               names_to = "sex",
               values_to = "death") %>% 
  pivot_wider(names_from = status,
              values_from = death) %>%
  replace(is.na(.), 0) %>% 
  mutate(Total = Yes + No)   

pop <- expand_grid(pop,
                   cause = unique(deaths$cause))

prevalence <- dt %>%
  dplyr::select(year, wt, id, age, sex:epilepsyev) %>%
  pivot_longer(-c(id:sex, year, wt),
               names_to = "cause",
               values_to = "status") %>%
  filter(status == "Yes") %>%
  group_by(year, sex, cause, age) %>% 
  summarise(prev = sum(wt), .groups ="drop") %>% 
  pivot_wider(names_from = sex,
              values_from = prev) %>% 
  mutate(Total = Male + Female) %>% 
  pivot_longer(c(Male:Total),
               names_to = "sex",
               values_to = "prev")


result <- function(.data) {
  regbsp <- gam(prev ~ s(age, year) + offset(log(pop)),
                data = .data,
                family = quasipoisson)
  
  model <-  function(a, y) {
    
    predict(regbsp, newdata = data.frame(age = a,
                                         year = y,
                                         pop = 1))
    
  }
  
  final <- outer(1:85, 2010:2018, model) %>% 
    exp() %>% 
    as.data.frame() %>% 
    mutate(age = 1:85) %>% 
    set_names(c(2010:2018, "age")) %>% 
    pivot_longer(-age,
                 names_to = "year",
                 values_to = "val") %>% 
    mutate(year = as.numeric(year))
  
  return(final)
  
}

full_dt_prev <- pop %>%
  full_join(prevalence) %>%
  replace(is.na(.), 0) %>%
  group_nest(sex, cause) %>% 
  mutate(model = map(data, ~ .x %>% 
                       result())) %>% 
  mutate(data = map(data, ~ .x %>% 
                      mutate(prv = prev / pop) %>% 
                      dplyr::select(-prev, -pop))) %>% 
  mutate(prg = map2(.x = model, .y = data, ~ .y %>%
                      full_join(.x))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(prg)

pop1 <- full_dt_prev %>%
  full_join(pop) %>%
  mutate(pop_u = val * pop,
         pop_h = pop - pop_u) %>% 
  dplyr::select(-c(prv, val))



# %>% 
#   filter(year == 2014) %>%
#   ggplot() +
#   # geom_line(aes(x = age, y = pop)) +
#   geom_line(aes(x = age, y = pop_u), linetype = 2) +
#   facet_wrap(cause ~ sex, scales = "free_y") + 
#   theme_bw()


# ok
full_dt_prev %>% 
  filter(year == 2014) %>%
  ggplot() +
  geom_line(aes(x = age, y = val, color = sex)) +
  geom_line(aes(x = age, y = prv, color = sex), linetype = 2) +
  facet_wrap( ~ cause, scales = "free_y") + 
  theme_bw() + 
  theme(legend.position = "bottom")

full_dt_mort <- pop1 %>%
  full_join(deaths) %>%
  replace(is.na(.), 0) %>%
  mutate(mx = Total / pop,
         mx_u = Yes / pop_u,
         mx_h = No / pop_h) %>% 
  dplyr::select(-c(contains("pop"), Yes, No, Total)) %>%
  pivot_longer(-c(sex:age),
               names_to = "status",
               values_to = "mx") %>%
  group_nest(sex, cause, year, status) %>%
  mutate(smt = map(data, ~ loess(.x$mx ~ .x$age)$fitted %>% 
                     as.data.frame() %>%
                     set_names("val") %>% 
                     mutate(age = 1:85,
                            val = ifelse(val < 0, 0 , val)))) %>%
  mutate(prg = map2(.x = smt, .y = data, ~ .y %>%
                      full_join(.x))) %>% 
  dplyr::select(-data, -smt) %>% 
  unnest(prg)



# ok
full_dt_mort %>%
  filter(year == 2012) %>%
  filter(status == "mx") %>% 
  ggplot() +
  geom_line(aes(x = age, y = mx, color = sex), alpha = .3) +
  geom_line(aes(x = age, y = val, color = sex), linetype = 2) +
  facet_wrap( ~ cause, scales = "free_y") + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  scale_y_log10()

# create final dataset
mrt <- full_dt_mort %>% 
  dplyr::select(-mx) %>% 
  pivot_wider(names_from = status,
              values_from = val) %>% 
  full_join(
    dplyr::select(full_dt_prev)
  ) %>% 
  dplyr::select(sex, cause, year, age, mx_u, mx_h, mx, prev = prv, prev_hat = val)
mrt |> write_csv("nhis_results.csv")
# save(mrt, file = "results.RData")