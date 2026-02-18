library(ipumsr)
library(tidyverse)
library(janitor)
library(mgcv)
library(MortalitySmooth)

ddi         <- read_ipums_ddi("nhis_00002.xml")
data        <- read_ipums_micro(ddi)
names(data) <- tolower(names(data))

fit_data <- expand_grid(age = 18:85,
                        pop = 1)

# HMD

hmd_ded <- HMDHFDplus::readHMDweb("USA","Deaths_1x1", 
                                  username = Sys.getenv("us"), 
                                  password = Sys.getenv("pw")) %>% 
  filter(between(Year, 2010, 2018)) %>% 
  dplyr::select(-Total) %>% 
  pivot_longer(c(Male, Female), 
               names_to = "sex",
               values_to = "dx") %>%
  set_names(tolower(names(.)))

hmd_pop <- HMDHFDplus::readHMDweb("USA","Exposures_1x1", 
                                  username = Sys.getenv("us"), 
                                  password = Sys.getenv("pw")) %>% 
  filter(between(Year, 2010, 2018)) %>% 
  dplyr::select(-Total) %>% 
  pivot_longer(c(Male, Female), 
               names_to = "sex",
               values_to = "pop") %>%
  set_names(tolower(names(.)))

hmd <- hmd_ded %>% 
  full_join(hmd_pop) %>% 
  mutate(age = ifelse(age > 84, 85, age)) %>% 
  group_by(year, sex, age) %>% 
  summarise(dx = sum(dx),
            pop = sum(pop),
            .groups = "drop")

# Smooth HMD mortality using P-splines approach of Camarda 2012.
# lambda = 100 is quite low, so we should preserve major features;
# any anomalies <20 or <90 aren't problematic for our analysis
hmd_s <-
  hmd |> 
  group_by(year, sex) |> 
  mutate(mx_raw = dx / pop,
         mx = Mort1Dsmooth(x = age, y = dx, lambda = 100,offset = log(pop), method = 3)$logmortality |> exp())


# problematic ratios are truncated out anyway for our analysis
hmd_s |> 
  mutate(ratio = mx / mx_raw) |> 
  ggplot(aes(x = age, y= ratio, color = as.factor(year))) + 
  geom_line() +
  facet_wrap(~sex) +
  scale_y_log10()

hmd_s %>% 
  filter(year %% 5 == 0) |> 
  ggplot(aes(x = age, y = mx))+
  geom_line()+ 
  facet_wrap(year~sex) +
  scale_y_log10() +
  geom_point(mapping = aes(y = mx_raw), alpha = .3)

# mortality
nhis_dat <- data %>%
  filter(mortelig == 1)

deaths <- nhis_dat %>%
  # calculate the age at death
  mutate(death_age = ifelse(mortstat == 1, 
                            mortdody - (year - age),
                            2019 - (year - age)), 
         # keep only people who died
         d.event = ifelse(mortstat == 1, 1, 0)) %>% 
  # create a new mortality weight
  mutate(m_wt = coalesce(mortwt, mortwtsa)) %>% 
  filter(mortdody < 2019) %>% 
  filter(d.event == 1)

dead <- deaths %>%
  dplyr::select(year,
                age,
                sex,
                cancerev:hyp2time,
                mortdody,
                death_age,
                d.event,
                m_wt) %>%
  dplyr::select(-c(dementiaev, epilepsyev)) %>%
  filter(d.event == 1) %>% 
  # remove people who do died 2 years after survey
  # respect the quarters too
  mutate(keep = ifelse(
    (mortdody - year < 2), 1, 0
  )) %>% 
  filter(keep == 1) %>%
  as_factor() %>%
  mutate(age = as.character(age)) %>%
  mutate(age = parse_number(age)) %>%
  mutate(mortdody = as.numeric(as.character(mortdody))) %>%
  pivot_longer(c(cancerev:hyp2time),
               names_to = "cause",
               values_to = "val") %>%
  dplyr::select(-c(d.event, keep)) %>% 
  filter(val %in% c("Yes", "No")) %>% 
  mutate(death_age = ifelse(death_age > 84, 85, death_age)) %>% 
  # remove NIU people
  group_by(mortdody, sex, death_age, cause, val) %>% 
  summarise(d = sum(m_wt), .groups = "drop") %>% 
  rename(age = death_age, year = mortdody)

pop1 <- data %>% 
  dplyr::select(year,
                age,
                sex,
                cancerev:hyp2time,
                perweight) %>% 
  filter(year < 2019) %>% 
  dplyr::select(-c(dementiaev, epilepsyev)) %>% 
  as_factor() %>%
  mutate(age = as.character(age)) %>%
  mutate(age = parse_number(age)) %>%
  pivot_longer(c(cancerev:hyp2time),
               names_to = "cause",
               values_to = "val") %>% 
  filter(val %in% c("Yes", "No")) %>% 
  # remove NIU people
  group_by(year, sex, age, cause, val) %>% 
  summarise(pop = sum(perweight), .groups = "drop")




# data for predict
dt_fit <- expand_grid(year = unique(pop1$year),
                      age = 18:85,
                      pop = 1)

dt_fit_pop <- expand_grid(year = unique(pop1$year),
                          age = 18:85)

# test this concept
dx <- hmd_s |> filter(year == 2015, sex == "Female") |> pull(dx)
pop <- hmd_s |> filter(year == 2015, sex == "Female") |> pull(pop)
age <- hmd_s |> filter(year == 2015, sex == "Female") |> pull(age)



mort1d_sm_pred <- function(age, dx, pop){
  
  dx[is.na(dx)] = 0

  mod = Mort1Dsmooth(x = age, 
                     y = dx, 
                     lambda = 5e4,
                     offset = log(pop), 
                     method = 3)
  mx  = predict(mod, age = 18:85) |> exp()
  # plot(18:85, mx, log = 'y')
  out = tibble(age = 18:85, mx = mx)
  out
}

# smoothing function
smoothit <- function(.data) { 
  
  model <- gam(d ~ s(age, bs = "ps") + year + offset(log(pop)),
               data = .data,
               family = quasipoisson)
  
  # Predict the smoothed population counts
  smoothed_pop <- predict(model_p, dt_fit_pop, type = "response")
  pred         <- predict(model,   dt_fit,     type = "response")
  
  final <- dt_fit %>%
    mutate(mx  = pred,
           pop = smoothed_pop)
  
  return(final)
  
}

pop1 |> 
  full_join(dead, by = join_by(year, sex, age, cause, val)) |> 
  #filter(sex == "Female", cause == "cancerev", year == 2015, val == "No")
  group_nest(sex, cause, val, year) |> 
  reframe(~mort1d_sm_pred(age = age,dx=d,pop=pop))

full <- pop1 %>% 
  full_join(dead) %>%
  group_nest(sex, cause, val) %>% 
  mutate(data = map(data, ~ .x %>% 
                      smoothit)) %>% 
  unnest(data)

full %>% 
  filter(sex == "Male", year == 2010) %>% 
  ggplot(aes(x = age, y = pop, color = cause)) + 
  geom_line() +
  scale_y_log10()+
  facet_wrap(~ val)


z <- full %>%
  full_join(pop1) %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  dplyr::select(-pop_No) %>% 
  rename(prev = pop_Yes) 
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(
    Ra         = mx_Yes / mx_No,
    mx_yes_new = mx / (1 - prev + prev * Ra),
    mx_no_new  = mx_No * Ra
  )



full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>% 
  full_join(hmd, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_yes_new = mx / (1 - prev + prev * Ra),
         mx_no_new  = mx_yes_new * Ra) %>% 
  filter(sex == "Male", year == 2010) %>% 
  ggplot() + 
  geom_line(aes(x = age, y = prev, color = cause)) 




full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>% 
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_yes_new = mx / (1 - prev + prev * Ra),
         mx_no_new  = mx_yes_new * Ra) %>% 
  filter(sex == "Male", year == 2010, cause == "cancerev") %>% 
  ggplot() + 
  geom_line(aes(x = age, y = mx_yes_new, color = cause)) +
  geom_line(aes(x = age, y = mx_no_new, color = cause), lty = 2) +
  scale_y_log10() + 
  #geom_line(aes(x = age, y = mx_Yes, color = cause), lty = 2)
  
  
  
  
  full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>%
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_no_new = mx / (1 - prev + prev * Ra),
         mx_yes_new  = mx_no_new * Ra) %>%
  filter(sex == "Male", year == 2010, cause == "cancerev") %>%
  pivot_longer(contains("mx"),
               names_to = "var",
               values_to = "val") %>%
  mutate(var = str_replace_all(var, pattern = "mx_", replacement = "mx-")) %>%
  filter(var != "mx") %>% 
  separate(var, c("one", "two"), sep = "-") %>%
  separate(two, c("two1", "three"), sep = "_") %>%
  mutate(three = ifelse(is.na(three), "old", three)) %>% 
  mutate(two1 = str_to_title(two1)) %>%
  filter(sex == "Male", year == 2010, cause == "cancerev") %>%
  ggplot(aes(x = age, y = val, color = two1, lty = three)) + 
  geom_line(linewidth = 1) +
  scale_y_log10()

heartattev  
cancerev
strokev

full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>%
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_no_new = mx / (1 - prev + prev * Ra),
         mx_yes_new  = mx_no_new * Ra) %>%
  filter(sex == "Male", year == 2010, cause == "heartattev") %>% 
  ggplot(aes(x = age, y = Ra, color = as.factor(year))) + 
  geom_line() + 
  facet_wrap(~ sex)


full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>%
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_no_new = mx / (1 - prev + prev * Ra),
         mx_yes_new  = mx_no_new * Ra) %>%
  filter(sex == "Male", year == 2010, cause == "strokev") %>%
  pivot_longer(contains("mx"),
               names_to = "var",
               values_to = "val") %>%
  mutate(var = str_replace_all(var, pattern = "mx_", replacement = "mx-")) %>%
  filter(var != "mx") %>% 
  separate(var, c("one", "two"), sep = "-") %>%
  separate(two, c("two1", "three"), sep = "_") %>%
  mutate(three = ifelse(is.na(three), "old", three)) %>% 
  mutate(two1 = str_to_title(two1)) %>%
  filter(sex == "Male", year == 2010, cause == "strokev") %>%
  ggplot(aes(x = age, y = val, color = two1, lty = three)) + 
  geom_line(linewidth = 1) +
  scale_y_log10()




full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  mutate(prev = pop_Yes / (pop_Yes + pop_No)) %>%
  full_join(hmd_s, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_no_new = mx / (1 - prev + prev * Ra),
         mx_yes_new  = mx_no_new * Ra) %>%
  filter(sex == "Male", year == 2010, cause == "hypertenev") %>%
  pivot_longer(contains("mx"),
               names_to = "var",
               values_to = "val") %>%
  mutate(var = str_replace_all(var, pattern = "mx_", replacement = "mx-")) %>%
  filter(var != "mx") %>% 
  separate(var, c("one", "two"), sep = "-") %>%
  separate(two, c("two1", "three"), sep = "_") %>%
  mutate(three = ifelse(is.na(three), "old", three)) %>% 
  mutate(two1 = str_to_title(two1)) %>%
  filter(sex == "Male", year == 2010, cause == "hypertenev") %>%
  ggplot(aes(x = age, y = val, color = two1, lty = three)) + 
  geom_line(linewidth = 1) +
  scale_y_log10()





full %>%
  pivot_wider(names_from = val,
              values_from = c(mx, pop)) %>%
  dplyr::select(-pop_No) %>% 
  rename(prev = pop_Yes) %>%
  full_join(hmd, by = c("year", "age", "sex")) %>% 
  mutate(Ra = mx_Yes / mx_No,
         mx_yes_new = mx / (1 - prev + prev * Ra),
         mx_no_new = mx_No * Ra) %>% 
  filter(sex == "Male", year == 2010, cause != "hyp2time") %>%
  pivot_longer(contains("mx"),
               names_to = "var",
               values_to = "val") %>% 
  ggplot(aes(x = age, y = val, color = var)) + 
  geom_line() +
  scale_y_log10() +
  facet_wrap(~ cause)





z %>% 
  filter(sex == "Male", year == 2010) %>% 
  ggplot() + 
  geom_line(aes(x = age, y = mx_No, color = cause)) +
  geom_line(aes(x = age, y = mx_Yes, color = cause), lty = 2) +
  scale_y_log10()






z %>% 
  