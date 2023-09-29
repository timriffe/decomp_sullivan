
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
                   control = list(lambda = 1e6, deg = 3,kr=8))$fitted * offi
  
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
  filter(age <= 102)

dec_data |> 
  ggplot(aes(x = age, y = prev, color = sex, linetype = country)) +
  geom_line() 

# these prevalence age patterns seem to reach values that are too high,
dec_data <- 
dec_data |> 
  mutate(prev = if_else(prev < 0.0001,0,prev))

dec_data_expanded <- list()
age <- 0:102
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
  theme_minimal() +
  xlim(50,100) +
  scale_y_log10()

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

library(xtable)
dec_data_expanded |> 
  filter(between(age,40,100)) |> 
  group_by(country, sex) |> 
  summarize(HLE = sully_normal(mx, prev, type = 'h'),
            ULE = sully_normal(mx, prev, type = "u")) |> xtable()

dec_data_expanded |> 
  filter(between(age,40,100)) |> 
  group_by(country, sex) |> 
  summarize(HLE = sully_rates(qhx = qhx, 
                              qux = qux, 
                              phux = phux, 
                              p0 = prev[1], 
                              type = 'h'),
            ULE = sully_rates(qhx = qhx, 
                              qux = qux, 
                              phux = phux, 
                              p0 = prev[1], 
                              type = 'u'))
# good match
write_csv(dec_data_expanded,"PAA_abstract_dec_data_expanded.csv")
options(scipen = 5)
fig1 <-
dec_data_expanded |> 
  # mutate(R_guess = rep(R_guess,4)) |> 
  filter(age<=100) |> 
  select(country, sex, age, `\\pi(a)` = prev, `m(a)` = mx) |> 
  pivot_longer(c(`\\pi(a)`,`m(a)`), names_to = "input", values_to = "value") |> 
  ggplot(aes(x=age,y = value, color = input, linetype = sex)) +
  geom_line() +
  scale_y_log10() +
  xlim(40,100) +
  facet_wrap(~country) +
  theme_minimal() +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16))  +
  labs(y = "log(value)")
  ggsave(fig1,file="fig1.pdf")

fig2 <- 
  tibble(age = age, `R(a)` = R_guess) |> 
  ggplot(aes(x=age, y = `R(a)`)) +
  geom_line() +
    theme_minimal()+theme(axis.title = element_text(size=16),
                          axis.text = element_text(size=14),
                          strip.text = element_text(size=16),
                          legend.text = element_text(size = 14),
                          legend.title = element_text(size=16)) +
    xlim(40,100) +
    ylim(0,8)
ggsave(fig2, file="fig2.pdf")

fig3 <-
  dec_data_expanded |> 
  filter(age <= 100) |> 
  mutate(`q(a)` = mx_to_qx(mx)) |> 
  select(country, sex, age, `qu(a)` = qux, `qh(a)` = qhx,`q(a)`,`phu(a)` = phux  ) |> 
  pivot_longer(-c(country, sex,age), values_to = "value", names_to = "component") |> 
  ggplot(aes(x=age, y = value, color = component, linetype = sex)) +
  geom_line() +
  facet_wrap(~country) +
  theme_minimal() +
  scale_y_log10() +
  xlim(40,100)+
  theme_minimal()+
  labs(y="transition") +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16))
fig3
ggsave(fig3, file = "fig3.pdf")

ctry_list <- list()
for (ctry in c("Germany", "United States of America")) {
  m <- dec_data_expanded |>
    filter(country == ctry, sex == "Male",
           between(age,40,100))
  f <- dec_data_expanded |>
    filter(country == ctry, sex == "Female",
           between(age,40,100))
  
  type_list <- list()
  for (type in c("h", "u", "t")) {
    method_list <- list()
    for (method in c("Sullivan", "Incidence")) {
      if (method == "Incidence") {
        cci <- sully_rates_decomp(
          mx1 = m$mx,
          mx2 = f$mx,
          pux1 = m$prev,
          pux2 = f$prev,
          R1 = R_guess[(41:101)],
          R2 = R_guess[(41:101)],
          # R1 = 5.9,
          # R2 = 5.9,
          type = type
        )  |>
          mutate(method = method) |>
          mutate(expectancy = type, .before = 1) |>
          mutate(country = ctry, .before = 1)
      }
      
      if (method == "Sullivan") {
        cci <- sully_normal_decomp(
          mx1 = m$mx,
          mx2 = f$mx,
          pux1 = m$prev,
          pux2 = f$prev,
          type = type
        )  |>
          mutate(method = method) |>
          mutate(expectancy = type, .before = 1) |>
          mutate(country = ctry, .before = 1)
        
      }
      method_list[[method]] <- cci
    }
    type_list[[type]] <- bind_rows(method_list)
  }
  ctry_list[[ctry]] <- bind_rows(type_list)
}
dec_results <- bind_rows(ctry_list)

components <-
dec_results |> 
  group_by(country, method, expectancy, component) |> 
  summarize(contribution  = sum(value)) 
components |> 
  pivot_wider(names_from = c(method, component), names_sep = " ", 
              values_from = contribution) |> 
  xtable()

fig4 <- 
components |> 
  filter(method == "Sullivan") |> 
  mutate(component = case_when(component == "mx" ~ "m(a)",
                               component == "pux" ~ "\\pi(a)")) |> 
  ggplot(aes(x = component, y = contribution, fill = component))+
  geom_col() +
  facet_grid(vars(expectancy),vars(country)) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16)) +
  guides(fill = "none") +
  scale_fill_discrete_qualitative("Warm")
ggsave(fig4, file = "fig4.svg")


fig5 <- 
components |> 
  filter(method == "Incidence") |> 
  mutate(component = case_when(component == "init" ~ "\\pi(40)",
                               component == "qhx" ~ "qh(a)",
                               component == "qux" ~ "qu(a)",
                               component == "onset" ~ "phu(a)")) |> 
  ggplot(aes(x = component, y = contribution, fill = component))+
  geom_col() +
  facet_grid(vars(expectancy),vars(country)) +
  theme_minimal() +
  geom_hline(yintercept = 0)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=16)) +
  guides(fill = "none") +
  scale_fill_discrete_qualitative("Dark3")
ggsave(fig5,file = "fig5.svg", width=8,height=7,units = "in")
library(colorspace)
hcl_palettes(plot=TRUE)
