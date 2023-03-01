# X <- model.matrix(formula_m_bart, data = meps)
# Y <- scale(meps$logY)
# M <- scale(meps$phealth)
# A <- meps$smoke
# seed <- 111111
rlasso_mediate <- function(X, Y, M, A, seed) {
  library(rlearner)
  
  set.seed(seed)
  rlearn_m <- rlasso(x = X, w = A, y = M, k_folds = 30)
  
  set.seed(seed)
  rlearn_total <- rlasso(
    x = X,
    y = Y,
    w = A,
    k_folds = 30
  )
  
  set.seed(seed)
  rlearn_indirect <- rlasso(
    x = cbind(X, A),
    y = Y,
    w = M,
    k_folds = 30
  )
  
  total <- rlearn_total$tau_hat %>% as.numeric()
  tau <- rlearn_m$tau_hat %>% as.numeric()
  delta <- tau * (rlearn_indirect$tau_hat %>% as.numeric())
  d <- rlearn_indirect$tau_hat %>% as.numeric()
  zeta <- total - delta
  
  return(list(total = total, delta = delta, tau = tau, d = d, zeta = zeta))
  
}
