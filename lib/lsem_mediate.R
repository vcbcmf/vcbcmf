lsem_mediate <- function(fit_m, fit_y, formula_X,
                         data_test, mediator_name, treat_name) {

  X <- model.matrix(formula_X, data = data_test)

  gamma_m <- coef(fit_m)[str_detect(names(coef(fit_m)), treat_name)]
  gamma_y <- coef(fit_y)[str_detect(names(coef(fit_y)), treat_name)]
  xi_y <- coef(fit_y)[str_detect(names(coef(fit_y)), mediator_name)]

  tau_hat <- as.numeric(X %*% gamma_m)
  zeta_hat <- as.numeric(X %*% gamma_y)
  d_hat <- as.numeric(X %*% xi_y)
  delta_hat <- tau_hat * d_hat

  sigma_m <- sigma(fit_m)
  sigma_y <- sigma(fit_y)

  return(list(zeta = zeta_hat,
              delta = delta_hat,
              tau_hat = tau_hat,
              sigma_m = sigma_m,
              sigma_y = sigma_y))
}