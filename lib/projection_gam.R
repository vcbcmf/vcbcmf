projection_gam <- function(data, samples) {
  projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
  for (i in 1:nrow(samples)) {
    formula <- samples[i,] ~ s(age, k = 10) + s(bmi, k = 10) + edu +
      s(log(income + 1000), k = 10) + s(povlev, k = 10) + region + sex + 
      marital + race + seatbelt
    gam_fit <- gam(formula, data = data)
    projections[i,] <- predict(gam_fit, type = 'response')
  }
  return(projections)
}