simulate_rlearn_med <- function(fit_m, fit_tau, fit_delta, w, sigma_y, sigma_m)
{
  
  ## Mediator Quants
  M_hat <- fit_m$m_hat %>% as.numeric()
  tau_m_hat <- fit_m$tau_hat %>% as.numeric()
  p_m_hat <- fit_m$p_hat %>% as.numeric()
  
  ## Mediator Draw
  M <- M_hat + (w - p_m_hat) * tau_m_hat + sigma_m * rnorm(length(M_hat))
  
  ## Outcome Quants
  Y_hat <- fit_tau$m_hat %>% as.numeric()
  A_hat <- fit_tau$p_hat
  tau_y <- fit_tau$tau_hat %>% as.numeric()
  d <- fit_delta$tau_hat %>% as.numeric()
  delta <- d * tau_m_hat
  zeta <- tau_y - delta
  
  Y <- Y_hat +
    (w - A_hat) * zeta +
    (M - M_hat) * d +
    sigma_y * rnorm(length(M_hat))
  
  return(list(Y = Y, M = M, delta = delta, zeta = zeta, tau = tau_y))
  
}