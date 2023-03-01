lsem_boot <- function(fit_m, fit_y, formula_X, data_train, data_test,
                      mediator_name, treat_name, outcome_name) {
  
  formula_m <- formula(fit_m)
  formula_y <- formula(fit_y)
  
  m_fitted <- predict(fit_m)
  m_resid <- sample(resid(fit_m), replace = TRUE)
  m_boot <- m_fitted + m_resid
  data_train_boot_m <- data_train
  data_train_boot_m[[mediator_name]] <- m_boot
  fit_m_boot <- lm(formula_m, data = data_train_boot_m)
  
  y_fitted <- predict(fit_y)
  y_resid <- sample(resid(fit_y), replace = TRUE)
  y_boot <- y_fitted + y_resid
  data_train_boot_y <- data_train
  data_train_boot_y[[outcome_name]] <- y_boot
  fit_y_boot <- lm(formula_y, data = data_train_boot_y)
  
  out_boot <- lsem_mediate(fit_m_boot, fit_y_boot, formula_X,
                           data_test, mediator_name, treat_name)
  
  return(list(zeta = out_boot$zeta,
              delta = out_boot$delta))
}