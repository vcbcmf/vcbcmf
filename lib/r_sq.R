r_sq <- function(samples, samples_proj) {
  r2 <- rep(NA, nrow(samples))
  for (i in 1:nrow(samples)) {
    num <- sum((samples[i,] - samples_proj[i,])^2)
    denom <- sum((samples[i,] - mean(samples[i,]))^2)
    r2[i] <- 1 - (num / denom)
  }
  
  samples_mean <- colMeans(samples)
  samples_proj_mean <- colMeans(samples_proj)
  num <- sum((samples_mean - samples_proj_mean)^2)
  denom <- sum((samples_mean - mean(samples))^2)
  r2_mean <- 1 - (num / denom)
  
  return(list(r2 = r2, r2_mean = r2_mean))
}

r_sq_gam_tree <- function(data, formula_y, samples) {
  samples_mean <- colMeans(samples)
  X <- model.matrix(formula_y, data = data)
  formula <- samples_mean ~ .
  
  ## Rsq for tree
  tree_fit <- rpart(formula, data = as.data.frame(X))
  projections <- predict(tree_fit)
  num <- sum((samples_mean - projections)^2)
  denom <- sum((samples_mean - mean(samples_mean))^2)
  r_sq_tree <- 1 - num / denom
  
  ## Rsq for Gam
  formula <- samples_mean ~ s(age, k = 10) + s(bmi, k = 10) + edu +
    s(log(income + 1000), k = 10) + s(povlev, k = 10) + region + sex +
    marital + race + seatbelt
  gam_fit <- gam(formula, data = data)
  projections <- predict(gam_fit, type = 'response')
  num <- sum((samples_mean - projections)^2)
  denom <- sum((samples_mean - mean(samples_mean))^2)
  r_sq_gam <- 1 - num / denom
  
  
  return(return(list(r_sq_tree = r_sq_tree, r_sq_gam = r_sq_gam)))
}
