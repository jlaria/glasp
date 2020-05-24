new_linear_classification <- function(beta, intercept, clusters, info, submodel, blueprint) {
  hardhat::new_model(beta = beta,
                     intercept = intercept,
                     clusters = clusters,
                     info = info,
                     blueprint = blueprint,
                     submodel = submodel,
                     class = "linear_classification")
}
