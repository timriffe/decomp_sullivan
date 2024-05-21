# install.packages("ipumsr")
library(ipumsr)
library(tidyverse)
library(janitor)
library(mgcv)

# R.utils::gunzip("nhis_00002.dat.gz", remove = FALSE)
ddi         <- read_ipums_ddi("dat_ipums/nhis_00002.xml")
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
                mortstat, mortdodq, mortdody, 
                wt_m1 = mortwt, wt_m2 = mortwtsa, 
                sampweight) %>% 
  as_factor() %>% 
  filter(year < 2019) %>% 
  dplyr::select(-c(dementiaev, epilepsyev))

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
  dplyr::select(-c(mortdodq)) %>%
  # remove missing data
  filter(!age %in% c("Unknown-don't know", "Unknown-refused")) %>% 
  mutate(age = ifelse(str_detect(age, "85"), "85", age),
         age = as.numeric(age)) %>%
  # this period has almost no information on dementia
  # remove year with missing weights
  mutate(year = as.numeric(as.character(year))) %>% 
  # mortality weight
  mutate(m_wt = coalesce(wt_m1, wt_m2))

# overall population
pop <- dt %>%
  dplyr::select(year, id, age, sex, wt) %>%
  # pivot_longer(c(cancerev:hyp2time),
  #              names_to  = "cause",
  #              values_to = "status") %>%
  # filter(status %in% c("Yes", "No")) %>%
  group_by(year, sex, age) %>%
  summarise(pop = sum(wt), .groups = "drop") %>%
  group_by(year, age) %>%
  # slow working but code looks nice
  # I can think of at least 2 identical ways of doing it faster
  # but the code will be longer
  group_modify(~ .x %>%
                 janitor::adorn_totals("row", name = "Total")) %>%
  ungroup()
# pivot_wider(names_from  = cause, 
#             values_from = pop,
#             values_fill = 0) 

# now lets calculate deaths
deaths <- dt %>%
  filter(mortstat == "Assumed deceased") %>%
  # year of death column is year now here
  dplyr::select(m_wt, id, age, sex:hyp2time, year = mortdody) %>%
  pivot_longer(-c(m_wt:sex, year),
               names_to = "cause",
               values_to = "status") %>%
  filter(status %in% c("No", "Yes")) %>% 
  group_by(year, sex, cause, status, age) %>% 
  summarise(death = sum(m_wt), .groups = "drop") %>%
  mutate(year = as.numeric(as.character(year))) %>%
  group_by(year, age, cause, status) %>%
  group_modify(~ .x %>%
                 janitor::adorn_totals("row", name = "Total")) %>% 
  ungroup() %>% 
  group_by(year, age, cause, sex) %>%
  group_modify(~ .x %>%
                 janitor::adorn_totals("row", name = "Total")) %>%
  mutate(sex    = as.character(sex),
         status = as.character(status))  

# expand population
pop <- expand_grid(pop,
                   cause = unique(deaths$cause)) %>% 
  mutate(sex = as.character(sex))



# calculate prevalence
prevalence <- dt %>%
  dplyr::select(year, wt, id, age, sex:hyp2time) %>%
  pivot_longer(-c(id:sex, year, wt),
               names_to = "cause",
               values_to = "status") %>%
  filter(status == "Yes") %>%
  group_by(year, sex, cause, age) %>% 
  summarise(prev = sum(wt), .groups ="drop") %>% 
  group_by(year, cause, age) %>% 
  group_modify(~ .x %>%
                 janitor::adorn_totals("row", name = "Total")) %>% 
  ungroup()

new_data <- expand_grid(age   = seq(1, 85, 1),
                        cause = unique(prevalence$cause),
                        year  = 2010:2018,
                        sex   = c("Male", "Female", "Total"),
                        pop   = 1)

full_dt_prev <- pop %>%
  full_join(prevalence, by = join_by(year, age, sex, cause)) %>%
  replace(is.na(.), 0) %>% 
  # filter(!is.na(prev)) %>%
  group_nest(sex, cause) %>% 
  nest_join(new_data, by = c("cause", "sex")) %>% 
  mutate(model = map(data, ~ gam(prev ~ s(age, bs = "ps") + 
                                   s(year, bs = "ps", k = 9) + offset(log(pop)),
                                 data = .,
                                 family = quasipoisson))) %>%
  mutate(pred = map2(.x = model, 
                     .y = new_data, ~ .y %>%
                       mutate(new = predict(.x, newdata = .y, type = "response"),
                              new = as.numeric(new)))) %>%
  mutate(data = map(data, ~ .x %>% 
                      mutate(prv = prev / pop) %>% 
                      dplyr::select(-prev, -pop))) %>% 
  mutate(prg = map2(.x = pred, .y = data, ~ .y %>%
                      full_join(.x, by = join_by(year, age)))) %>%
  dplyr::select(-c(data:pred)) %>% 
  unnest(prg) %>% 
  dplyr::select(-pop)

# check is ok
# most likely a data wrangling problem with this particular cause
# I remove it for now, can fix later
full_dt_prev %>%
  filter(year == 2018) %>%
  ggplot() +
  geom_line(aes(x = age, y = new, color = sex)) +
  geom_line(aes(x = age, y = prv, color = sex), linetype = 2) +
  facet_wrap( ~ cause, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom")

# healthy and unhealthy population for denominator for mx
# I calculate it following the Sullivan idea with the fitted prevalence
pop1 <- full_dt_prev %>%
  full_join(pop, by = join_by(sex, cause, year, age)) %>%
  mutate(Yes = new * pop,
         No = pop - Yes) %>%
  rename(Total = pop) %>% 
  dplyr::select(-c(new, prv)) %>%
  pivot_longer(c(Total:No),
               names_to  = "status",
               values_to = "pop")



# now for mortality smoothing
# were I had to go with the Local Polynomial Regression Fitting
# there is so many zeros in the data other models go crazy
# not a fast fix u unfortunately
new_data <- expand_grid(new_data,
                        status = c("Yes", "No", "Total"))


# 36 for more df than data
full_dt_mort <- pop1 %>%
  full_join(deaths, by = join_by(sex, cause, year, age,
                                 status)) %>%
  filter(year < 2019) %>%
  replace(is.na(.), 0) %>% 
  # filter(!is.na(death)) %>% 
  group_nest(sex, cause, status) %>%
  nest_join(new_data, by = c("cause", "sex", "status")) %>% 
  mutate(model = map(data, ~ gam(death ~ s(age, bs = "ps") + 
                                   s(year, bs = "ps", k = 9) + offset(log(pop)),
                                 data = .,
                                 family = quasipoisson))) %>%
  mutate(pred = map2(.x = model, 
                     .y = new_data, ~ .y %>%
                       mutate(new = predict(.x, newdata = .y, type = "response"),
                              new = as.numeric(new)))) %>%
  mutate(data = map(data, ~ .x %>% 
                      mutate(mx = death / pop) %>% 
                      dplyr::select(-death, -pop))) %>% 
  mutate(prg = map2(.x = pred, .y = data, ~ .y %>%
                      full_join(.x, by = join_by(year, age)))) %>%
  dplyr::select(-c(data:pred)) %>% 
  unnest(prg) %>% 
  dplyr::select(-pop)

full_dt_mort %>%
  filter(year == 2012) %>%
  filter(status == "Yes") %>%
  # filter(age > 30) %>%
  ggplot() +
  geom_line(aes(x     = age, 
                y     = new, 
                color = sex)) +
  geom_line(aes(x     = age, 
                y     = mx, 
                color = sex), 
            linetype  = 2) +
  facet_wrap( ~ cause, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom")

# create final dataset

full_dt <- full_dt_mort %>% 
  dplyr::select(-mx) %>% 
  pivot_wider(names_from  = status,
              values_from = new) %>% 
  full_join(
    dplyr::select(full_dt_prev, -prv)
  ) %>% 
  dplyr::select(sex, 
                cause,
                year,
                age,
                mx_u = Yes, 
                mx_h = No, 
                mx   = Total, 
                prev = new)

# save
full_dt|> 
  write_csv("nhis_results.csv")
