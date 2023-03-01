get_clever_cov <- function(data_train, data_test, formula_m,
                           mediator_name, outcome_name, treat_name) {
  
  X_m_test  <- model.frame(formula_m, data = data_test) %>%
    select(-all_of(mediator_name)) %>% 
    preprocess_df() %>% 
    pluck("X")
  
  X_m_train  <- model.frame(formula_m, data = data_train) %>%
    select(-all_of(mediator_name)) %>% 
    preprocess_df() %>% 
    pluck("X")
  
  X_m0_train <- X_m_train[data_train[[treat_name]] == 0,]
  X_m1_train <- X_m_train[data_train[[treat_name]] == 1,]
  m0_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 0]
  m1_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 1]
  
  bart_m0 <- softbart(X_m0_train, m0_train, X_m_test)
  bart_m1 <- softbart(X_m1_train, m1_train, X_m_test)
  
  m0_hat <- bart_m0$y_hat_test_mean
  m1_hat <- bart_m1$y_hat_test_mean
  
  return(list(m0_hat = m0_hat, m1_hat = m1_hat))
}
