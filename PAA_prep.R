library(tidyverse)
library(scales)
source("functions.R")
# Fig 1:




prev1 <- seq(0.01,.7,length=111)
prev2 <- seq(.01,.6,length=111)
a <- 0:110

mx1 <- .0002 * exp(.075 * a)
plot(mx1,log='y')
mx2 <- .0001 * exp(.075 * a)
plot(mx1, log = 'y')
lines(mx2)

mx<-
bind_rows(tibble(age = a, pop = "1", mx = mx1),
          tibble(age = a, pop = "2", mx = mx2))

prev <-
  bind_rows(tibble(age = a, pop = "1", prev = prev1),
            tibble(age = a, pop = "2", prev = prev2))

prev |> 
  ggplot(aes(x=age, y = prev, color = pop)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(y = "prevalence")

mx |> 
  ggplot(aes(x=age, y = mx, color = pop)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  scale_y_log10(labels = label_comma()) +
  labs(y = "mortality rate (log)")

mx |> 
  group_by(pop) |> 
  mutate(lx = mx_to_lx(mx) %>% '['(1:n())) |> 
  ggplot(aes(x=age, y = lx, color = pop)) +
  geom_line() +
  theme_minimal(base_size = 14)

f1 <-
mx |> 
  group_by(pop) |> 
  mutate(lx = mx_to_lx(mx) %>% '['(1:n())) |> 
  ggplot(aes(x=age, y = lx, color = pop)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  geom_line(data = prev, mapping = aes(y=prev), linetype = 2) +
  labs(y = "sullivan components") +
  annotate("text",40,.35, label = "prevalence",angle = 32,size=6) +
  annotate("text",75,.8, label = "survivorship",angle = -60, size = 6)


mx |> 
  group_by(pop) |> 
  mutate(lx = mx_to_lx(mx) %>% '['(1:n())) |> 
  rename(edad = age,
         pob = pop) |> 
  ggplot(aes(x=edad, y = lx, color = pob)) +
  geom_line(linewidth=1.5) +
  theme_minimal(base_size = 14) +
  geom_line(data = prev |> rename(edad=age, pob=pop), mapping = aes(y=prev), linetype = 2, linewidth=1.5) +
  labs(y = "Componentes Sullivan") +
  annotate("text",40,.35, label = "prevalencia",angle = 32,size=6) +
  annotate("text",75,.8, label = "supervivencia",angle = -60, size = 6)

f2 <-
mx |> 
  group_by(pop) |> 
  mutate(lx = mx_to_lx(mx) %>% '['(1:n())) |> 
  left_join(prev,by=join_by(pop,age)) |> 
  mutate(lhx = lx * (1-prev)) |> 
  ggplot(aes(x=age, y = lx, color = pop)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  geom_line(mapping = aes(y=prev), linetype = 2) +
  geom_line(mapping = aes(y=lhx), linetype = 5) +
  labs(y = "sullivan components")+
  annotate("text",40,.35, label = "prevalence",angle = 32,size=6,color=gray(.5)) +
  annotate("text",75,.8, label = "survivorship",angle = -60, size = 6,color=gray(.5)) +
  annotate("text",50,.72, label = "healthy years",angle = -50, size = 6) 

ggsave("paa1.png",f1, width=10,height=10,units="cm")
ggsave("paa2.png",f2, width=10,height=10,units="cm")


sully_normal(mx1,prev1,type='h') |> round(2)
sully_normal(mx2,prev2,type='h') |> round(2)

sully_normal_decomp(mx1, mx2, pux1=prev1, pux2=prev2, type= 'h') |> 
  group_by(component) |> 
  summarize(cc=sum(value))

# Check age margin
sully_normal_decomp(mx1, mx2, pux1=prev1, pux2=prev2, type= 'h') |> 
  group_by(age) |> 
  summarize(cc=sum(value)) |> 
  ggplot(aes(x = age, y = cc)) +
  geom_line() +
  theme_minimal()

Ra <- 2

case1 <- sully_derive_rates(mx1, prev1, 2)|> 
  mutate(mux = qx_to_mx(qux),
         mhx = qx_to_mx(qhx)) |> 
  select(age, mux, mhx, phux) |> 
  pivot_longer(mux:phux, names_to = "variable", values_to = "value") |> 
  mutate(case = "1", .before =1)
case2 <- sully_derive_rates(mx1, prev1, 4)|> 
  mutate(mux = qx_to_mx(qux),
         mhx = qx_to_mx(qhx)) |> 
  select(age, mux, mhx, phux) |> 
  pivot_longer(mux:phux, names_to = "variable", values_to = "value") |> 
  mutate(case = "2", .before =1)
# 
# test1 <- rates_derive_sully(qux = case1$qux,
#                    qhx = case1$qhx,
#                    phux = case1$phux,
#                    p0=0,
#                    age=0:110)
# test2 <- rates_derive_sully(qux = case2$qux,
#                             qhx = case2$qhx,
#                             phux = case2$phux,
#                             p0=0,
#                             age=0:110)

# These are not yet calibrated, has to do w inconsistencies here and there.
f4 <-
bind_rows(case1, case2) |> 
  ggplot(aes(x = age, y = value, color = variable, linetype = case)) +
  geom_line() +
  theme_minimal(base_size = 14)+
  scale_y_log10()


f3 <- tibble(age = 0:110, mx = mx1, prev = prev1) |> 
  pivot_longer(-age, values_to = "value", names_to = "variable") |> 
  ggplot(aes(x = age, y = value, color = variable)) +
  geom_line() +
  theme_minimal(base_size = 14) +
  scale_y_log10()

ggsave("paa3.png", f3, width = 10, height = 8, units = "cm")
ggsave("paa4.png", f4, width = 10, height = 8, units = "cm")

sully_rates_decomp(mx1,mx2,pux1 = prev1, pux2 = prev2, R1= 2, R2 = 2) |> 
  group_by(component) |> 
  summarize(cc = sum(value))

sully_rates_decomp(mx1,mx2,pux1 = prev1, pux2 = prev2, R1= 4, R2 = 4) |> 
  group_by(component) |> 
  summarize(cc = sum(value))

sully_rates_decomp(mx1,mx2,pux1 = prev1, pux2 = prev2, R1= 1, R2 = 1) |> 
  group_by(component) |> 
  summarize(cc = sum(value))

