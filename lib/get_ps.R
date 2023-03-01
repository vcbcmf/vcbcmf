get_ps <- function(data_train, data_test, formula_ps) {
                   glm_logit <- glm(formula_ps, data = data_train, family = binomial)
                   pi_hat <- predict(glm_logit, newdata = data_test,type = 'response')
                   return (pi_hat)
}