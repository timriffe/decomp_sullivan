
# note, step 7 returns an object with prev and mx,
# which we want to compare againt, externally
# provided prev and mx

# step 8

# needs 4 pieces for 2 transient states.
#N^2 for N transient states wuth full reversibility
pars <- c(slope = .08, intercept = .0001)
min_fun <- function(pars, mx_external, prev_external, model_pars, scale = 50){
  
  .data <- as_tibble(pars)
  stationary_outputs <- step7(.data)
  
  external_inputs <- tibble(age,
                            mx_external = mx_external,
                            prev_external = prev_external)
  compare <- left_join(stationary_outputs, external_inputs, by = join_by(age))
  
  RES1 <- 
  compare |> 
    mutate(mx_resid = abs(log(mx) - log(mx_external)),
           # maybe logit?
           prev_resid = abs(prev - prev_external)) |> 
    summarize(RES = sum(mx_resid + prev_resid)) |> 
    pull(RES)
  
  RES2 <- sum(scale * abs(pars - model_pars))
  
  RES1 + RES2
}

step8 <- function(model_pars, mx_external, prev_external, scale = 50){
  
  # change this to be 4 pars
  standard_pars <- c(slope1 = model_pars$slope, 
                  intercept2 = model_pars$intercept)
  
  opt_pars <- optim(par = standard_pars, 
                    fn = min_fun,
                    model_pars = standard_pars,
                    scale = scale)$par
  
  # convert to tibble
  
  # convert to transitions
  
  
}




