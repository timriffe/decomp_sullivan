
source("functions.R")

# TR: note I'm keeping for now all ages 0-110, and graduating to these too,
# I stopped the exercise on seeing the age patterns of prevalence, 
# which are impossible.

# read data


# Pop german females
pop_ger <- read_table("pop_germ.txt", 
                  skip = 2, 
                  show_col_types = FALSE) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(age = Age, Female, Male) %>% 
  pivot_longer(c(Female,Male), names_to = "sex", values_to = "pop") |> 
  mutate(age = parse_number(age)) %>%
  # mutate(age = if_else(age > 99, 100, age)) %>% 
  group_by(sex, age) %>% 
  summarise(pop = sum(pop), .groups = "drop")|> 
    mutate(country = "Germany")

# Pop USA males 
pop_usa <- read_table("pop_usa.txt", 
                   skip = 2, 
                   show_col_types = FALSE) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(age = Age, Female, Male) %>% 
  pivot_longer(c(Female,Male),names_to = "sex", values_to = "pop") |> 
  mutate(age = parse_number(age)) %>%
  # mutate(age = if_else(age > 99, 100, age)) %>% 
  group_by(sex,age) %>% 
  summarise(pop = sum(pop), .groups = "drop") |> 
  mutate(country = "United States of America")

pop <- bind_rows(pop_ger, pop_usa)
# deaths german females from all causes
dx_ger <- read_table("germany_deaths.txt", 
                     skip = 2,
                     show_col_types = FALSE) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(age = Age, Female, Male) %>% 
  pivot_longer(c(Female,Male),names_to = "sex",values_to = "deaths") |> 
  mutate(age = parse_number(age)) %>%
  # mutate(age = if_else(age > 99, 100, age)) %>% 
  group_by(sex, age) %>% 
  summarise(deaths = sum(deaths), .groups = "drop") |> 
  mutate(country = "Germany")

# deaths USA males from all causes
dx_usa <- read_table("usa_deaths.txt", 
                     skip = 2,
                     show_col_types = FALSE) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(age = Age, Female, Male) %>% 
  pivot_longer(c(Female,Male),
               names_to = "sex", 
               values_to = "deaths") |> 
  mutate(age = parse_number(age)) %>%
  #mutate(age = ifelse(age > 99, 100, age)) %>% 
  group_by(sex, age) %>% 
  summarise(deaths = sum(deaths), .groups = "drop") |> 
  mutate(country = "United States of America")

deaths <-
  bind_rows(dx_ger, dx_usa)

mort <- left_join(pop, deaths,by = join_by(sex, age, country)) |> 
  relocate(country, .before = 1)

# prevalence counts for usa and germany males and female 
# alzheimer's disease and other dementias
alz_counts <-
  read_csv("overall.csv", show_col_types = FALSE) %>% 
  dplyr::select(year, contains("name"), val) %>%
  set_names(str_remove(names(.), "_name$")) %>%
  filter(year == 2018, 
         cause == "Alzheimer's disease and other dementias") %>% 
  dplyr::select(country = location, sex, age, val) %>% 
  mutate(age = str_remove(age, " years$")) %>%
  mutate(age = str_remove(age, " year$")) %>% 
  mutate(age = case_when(age == "<1" ~ 0,
                         age == "95+" ~ 95,
                         age == "1-4" ~ 1,
                         age == "5-9" ~ 5,
                         TRUE ~ suppressWarnings(as.integer(substr(age,1,2))))) |>     mutate(val = ifelse(val == 0, 0.000001, val)) |> 
  arrange(country, sex, age)
alz_prev<-
  read_csv("IHME-GBD_2019_DATA-5b3d4823-1/IHME-GBD_2019_DATA-5b3d4823-1.csv") |> 
  filter(metric_name == "Percent") |> 
  select(country = location_name,
         sex = sex_name,
         age = age_name,
         val) |>   
  mutate(age = str_remove(age, " years$")) %>%
  mutate(age = str_remove(age, " year$")) %>% 
  mutate(age = case_when(age == "<1" ~ 0,
                         age == "95+" ~ 95,
                         age == "1-4" ~ 1,
                         age == "5-9" ~ 5,
                         TRUE ~ suppressWarnings(as.integer(substr(age,1,2))))) |>     mutate(val = ifelse(val == 0, 0.000001, val)) |> 
  arrange(country, sex, age)
pop_gbd <- pop |> 
  mutate(age2 = age - age %% 5,
         age2 = if_else(between(age,1,4),1,age2),
         age2 = if_else(age>95,95,age2)) |> 
  group_by(country, sex, age2) |> 
  summarize(pop = sum(pop), .groups = "drop") |> 
  rename(age = age2)

alz_counts <- pop_gbd |> 
  left_join(alz_prev,
            by = join_by(country, sex, age)) |> 
  mutate(count = val * pop)

# sanity check, very high closeout values I must say
# 60% in 95+, really? I don't believe it. 
  alz_prev|> 
    ggplot(aes(x=age,y=val,color=sex, linetype=country))+
    geom_step()

groups <- alz_counts |> 
  select(country, sex) |> 
  distinct()
# graduate to single age intervals
countsL <- list()
for (i in 1:nrow(groups)){
  meta <- groups[i, ]
  countsi <- 
    meta |> 
    left_join(alz_counts,
              by = join_by(country, sex))
  offi <- 
    meta |> 
    left_join(mort,
              by = join_by(country, sex)) |> 
    pull(pop)
   counts <- pclm( countsi$age,  
                   countsi$count, 
                   offset = offi, 
                   nlast = 16,
                   control = list(lambda = 1e5, deg = 3))$fitted * offi
   chunki <- cross_join(tibble(count = counts, age = 0:110), meta)
   countsL[[i]] <- chunki
}
alz_counts1 <- bind_rows(countsL)


# combine results, calculate the share with disability
dec_data <- left_join(mort, 
                      alz_counts1,
                      by = join_by(country, sex, age)) |> 
  mutate(mx = deaths / pop,
         prev = count / pop) |> 
  filter(age <= 100)

dec_data |> 
  ggplot(aes(x = age, y = prev, color = sex, linetype = country)) +
  geom_line() 

# these prevalence age patterns seem to reach values that are too high,
dec_data <- 
dec_data |> 
  mutate(prev = if_else(prev < 0.0001,0,prev))

dec_data_expanded <- list()

R_guess <- 7.964 + .01411 * age -.0002017 * age ^ 2 + -.000005974 * age^3


for (i in 1:nrow(groups)){
  meta <- groups[i, ]
  chunki <- dec_data |> 
    right_join(meta, 
               by = join_by(country, sex))
  ratesi <- sully_derive_rates(mx_all = chunki$mx,
                               pux = chunki$prev,
                               R_guess = R_guess)
  dec_data_expanded[[i]] <- 
    left_join(chunki, ratesi,by = join_by(age))
}
dec_data_expanded <- bind_rows(dec_data_expanded)

# check everything
dec_data_expanded |> 
  select(-c(pop,deaths,count)) |> 
  pivot_longer(mx:phux,names_to = "measure",values_to="value") |> 
  ggplot(aes(x = age, y = value, color = measure)) +
  geom_line() +
  facet_wrap(sex~country) +
  theme_minimal()

# just check rates
dec_data_expanded |> 
  select(-c(pop,deaths,count)) |> 
  mutate(qx_all = mx_to_qx(mx)) |> 
  pivot_longer(c(mx,qux,qhx),names_to = "measure",values_to="value") |> 
  ggplot(aes(x = age, y = value, color = measure, linetype = sex)) +
  geom_line() +
  facet_wrap(~country) +
  theme_minimal() +
  scale_y_log10()


dec_data_expanded |> 
  group_by(country, sex) |> 
  summarize(HLE = sully_normal(mx, prev, type = 'h'),
            ULE = sully_normal(mx, prev, type = "u"))

dec_data_expanded |> 
  group_by(country, sex) |> 
  summarize(HLE = sully_rates(qhx = qhx, 
                              qux = qux, 
                              phux = phux, 
                              p0 = 0, 
                              type = 'h'),
            ULE = sully_rates(qhx = qhx, 
                              qux = qux, 
                              phux = phux, 
                              p0 = 0, 
                              type = 'u'))
# good match

dec_results <- list()
for (i in nrow(groups)){
  meta <- groups[i, ]
  chunki <- dec_data_expanded |> right_join(meta)
}

ccr <- sully_rates_decomp(mx_all,
                          mx_all2,
                          pux,
                          pux2,
                          R1 = 1,
                          R2 = 1,
                          type = "h") 

ccr |> 
  ggplot(aes(x=age,y=value,color=component)) +
  geom_line() +
  geom_point(data = filter(ccr,component == "init"))



ccs <- sully_normal_decomp(
  mx_all,
  mx_all2,
  pux,
  pux2,
  type = "u")


ccs |> 
  ggplot(aes(x=age,y=value,color=component)) +
  geom_line()
