
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
                   countsi$val, 
                   offset = offi, 
                   nlast = 16)$fitted * offi
   chunki <- cross_join(tibble(count = counts, age = 0:110), meta)
   countsL[[i]] <- chunki
}
alz_counts1 <- bind_rows(countsL)


# combine results, calculate the share with disability
dec_data <- left_join(mort, 
                      alz_counts1,
                      by = join_by(country, sex, age)) |> 
  mutate(mx = deaths / pop,
         prev = count / pop)

dec_data |> 
  ggplot(aes(x = age, y = prev, color = sex, linetype = country)) +
  geom_line()

# these prevalence age patterns seem to reach values that are too high

# now to models
# females
mx_all  <- females$mx
pux     <- res$f

# males
mx_all2 <- males$mx
pux2    <- res$m
age     <- res$age



# good
sully_normal(mx_all, pux)
p1 <- sully_derive_rates(mx_all, pux, 1) # exact correspondence only have to find the R value
sully_rates(p1$qhx, p1$qux, p1$phux, p0 = pux[1])

# good
sully_normal(mx_all2, pux2)
p2 <- sully_derive_rates(mx_all2, pux2, 1)
sully_rates(p2$qhx, p2$qux, p2$phux, p0 = pux[1])



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
