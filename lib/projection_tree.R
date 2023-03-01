projection_tree <- function(data, formula_y, samples) {
  projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
  for (i in 1:nrow(samples)) {
    X <- model.matrix(formula_y, data = data)
    formula <- samples[i,] ~ .
    tree_fit <- rpart(formula, data = as.data.frame(X))
    projections[i,] <- predict(tree_fit)
  }
  return (projections)
}
