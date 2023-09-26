
source("functions.R")

# read data
# make 99 OA
# Pop german females
pop <- read_table("pop_germ.txt", skip = 2) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(Age, Female) %>% 
  mutate(Age = parse_number(Age)) %>%
  mutate(Age = ifelse(Age > 98, 100, Age)) %>% 
  group_by(Age) %>% 
  summarise(Female = sum(Female), .groups = "drop")

# Pop USA males 
pop2 <- read_table("pop_usa.txt", skip = 2) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(Age, Male) %>% 
  mutate(Age = parse_number(Age)) %>%
  mutate(Age = ifelse(Age > 98, 100, Age)) %>% 
  group_by(Age) %>% 
  summarise(Male = sum(Male), .groups = "drop")

# deaths german females from all causes
dx_ger <- read_table("germany_deaths.txt", skip = 2) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(Age, Female) %>% 
  mutate(Age = parse_number(Age)) %>%
  mutate(Age = ifelse(Age > 98, 100, Age)) %>% 
  group_by(Age) %>% 
  summarise(Female = sum(Female), .groups = "drop")

# deaths USA males from all causes
dx_usa <- read_table("usa_deaths.txt", skip = 2) %>% 
  filter(Year == 2018) %>% 
  dplyr::select(Age, Male) %>% 
  mutate(Age = parse_number(Age)) %>%
  mutate(Age = ifelse(Age > 98, 100, Age)) %>% 
  group_by(Age) %>% 
  summarise(Male = sum(Male), .groups = "drop")


# prevalence counts germany females for breast cancer
y <- read_csv("overall.csv") %>% 
  dplyr::select(year, contains("name"), val) %>%
  set_names(str_remove(names(.), "_name$")) %>%
  filter(year == 2018, sex == "Female", cause == "Breast cancer",
         location == "Germany") %>% 
  dplyr::select(age, val) %>% 
  mutate(age = str_remove(age, " years$")) %>%
  mutate(age = str_remove(age, " year$")) %>%
  mutate(age = factor(age, levels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24",
                                      "25-29", "30-34", "35-39", "40-44", "45-49", 
                                      "50-54", "55-59", "60-64", "65-74", "70-74", 
                                      "75-79", "80-84", "85-89", "90-94", "95+"))) %>% 
  arrange(age) %>% 
  mutate(age = as.character(age)) %>% 
  mutate(age = ifelse(age == "<1", "0", age)) %>%
  separate(age, c("one", "two"), sep = "-") %>%
  dplyr::select(age = one, val) %>% 
  mutate(age = ifelse(age == "95+", "95", age)) %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate(val = ifelse(val == 0, 0.000001, val))# pclm requirement


# prevalence counts USA males for Alzheimer's disease and other dementias
y2 <- read_csv("overall.csv") %>%
  dplyr::select(year, contains("name"), val) %>% 
  set_names(str_remove(names(.), "_name$")) %>% 
  filter(year == 2018, sex == "Male", cause == "Alzheimer's disease and other dementias",
         location == "United States of America") %>% 
  dplyr::select(age, val) %>% 
  mutate(age = str_remove(age, " years$")) %>%
  mutate(age = str_remove(age, " year$")) %>%
  mutate(age = factor(age, levels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24",
                                      "25-29", "30-34", "35-39", "40-44", "45-49", 
                                      "50-54", "55-59", "60-64", "65-74", "70-74", 
                                      "75-79", "80-84", "85-89", "90-94", "95+"))) %>% 
  arrange(age) %>% 
  mutate(age = as.character(age)) %>% 
  mutate(age = ifelse(age == "<1", "0", age)) %>%
  separate(age, c("one", "two"), sep = "-") %>%
  dplyr::select(age = one, val) %>% 
  mutate(age = ifelse(age == "95+", "95", age)) %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate(val = ifelse(val == 0, 0.000001, val))



# graduate to single age intervals
z  <- pclm( y$age,  y$val, nlast = 5)$fitted
z1 <- pclm(y2$age, y2$val, nlast = 5)$fitted


# combine results, calculate the share with disability
res <-  tibble(f = z,
               m = z1,
               age = 0:99, 
               pop_f = pop$Female, 
               pop_m = pop2$Male) %>% 
  mutate(m = m / pop_m,
         f = f / pop_f)



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
