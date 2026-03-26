
# isocline_slope_test.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

test_isocline_slope <- function(age, lx, prev, Rx,
                                age_int = 1,
                                bounds = c(1e-12, 2),
                                turnover_K = Inf,
                                hu_eps_log = 1e-2,
                                tiny = 1e-14,
                                line_tol = 1e-12){
  
  n <- length(age)
  stopifnot(length(lx)==n, length(prev)==n, length(Rx)==n)
  
  # stocks
  lH <- lx * (1 - prev)
  lU <- lx * prev
  
  # no‑returns hazards
  hz_nr <- derive_noreturns_hazards(age=age, lx=lx, prev=prev, Rx=Rx, age_int=age_int, verbose=FALSE)
  
  hd <- hz_nr |> dplyr::filter(trans=="hd") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  ud <- hz_nr |> dplyr::filter(trans=="ud") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  hu_nr <- hz_nr |> dplyr::filter(trans=="hu") |> dplyr::arrange(age) |> dplyr::pull(hazard)
  
  states <- build_states_targets(age, lx, prev)
  
  by_age <- lapply(seq_len(n), function(i){
    
    hu1 <- max(hu_nr[i], bounds[1])
    uh1 <- 0
    
    cap <- if(is.finite(turnover_K)) turnover_K * (hu1 + uh1) else Inf
    
    hu2 <- min(bounds[2], hu1 * exp(hu_eps_log))
    
    uh_hi <- min(bounds[2], if(is.finite(cap)) max(0, cap - hu2) else bounds[2])
    
    sol <- solve_uh_given_hu(i, states, hd, ud,
                             hu = hu2,
                             uh_lo = 0,
                             uh_hi = uh_hi,
                             age_int = age_int,
                             tiny = tiny)
    
    if(!is.finite(sol$loss) || sol$loss > line_tol){
      return(tibble::tibble(i=i, age=age[i], ok=FALSE,
                            m_pred=lH[i]/lU[i], m_emp=NA_real_,
                            rel_err=NA_real_, loss2=sol$loss))
    }
    
    uh2 <- sol$uh
    
    m_emp <- (uh2 - uh1)/(hu2 - hu1)
    m_pred <- lH[i]/lU[i]
    rel_err <- (m_emp - m_pred)/m_pred
    
    tibble::tibble(i=i, age=age[i], ok=TRUE,
                   m_pred=m_pred, m_emp=m_emp,
                   rel_err=rel_err, loss2=sol$loss)
  }) |> dplyr::bind_rows()
  
  overall <- by_age |> dplyr::summarise(
    n_ages=dplyr::n(),
    n_ok=sum(ok),
    median_abs_rel_err=median(abs(rel_err[ok]), na.rm=TRUE),
    max_abs_rel_err=max(abs(rel_err[ok]), na.rm=TRUE)
  )
  
  list(overall=overall, by_age=by_age)
}