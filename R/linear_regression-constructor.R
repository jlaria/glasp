new_linear_regression <- function(beta, intercept, clusters, info, blueprint) {
  hardhat::new_model(beta = beta,
                     intercept = intercept,
                     clusters = clusters,
                     info = info,
                     blueprint = blueprint,
                     class = "linear_regression")
}
