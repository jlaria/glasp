#' Predict from a `linear_classification`
#'
#' @param object A `linear_classification` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"prob"` for numeric predictions.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @export
predict.linear_classification <- function(object, new_data, type = "prob", ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_predict_types_classification())
  predict_linear_classification_bridge(type, object, forged$predictors)
}

valid_predict_types_classification <- function() {
  c("prob", "numeric")
}

# ------------------------------------------------------------------------------
# Bridge

predict_linear_classification_bridge <- function(type, model, predictors) {
  #predictors <- as.matrix(predictors)

  # Check submodel
  subtype <- paste0(type, "_", model$submodel)

  predict_function <- get_predict_function_classification(subtype)
  predictions <- predict_function(model, predictors)

  hardhat::validate_prediction_size(predictions, predictors)

  predictions
}

get_predict_function_classification <- function(type) {
  switch(
    type,
    prob_survival = predict_linear_classification_prob_survival,
    numeric_survival = predict_linear_classification_numeric_survival
  )
}

# ------------------------------------------------------------------------------
# Implementation

#' @importFrom "stats" "approx"
predict_linear_classification_prob_survival <- function(model, predictors) {
  # 1) We need the baseline survival function model$blin_S
  # 2) predictors has to be in the expanded format
  time <- predictors$time

  predictors$time <- NULL

  X <- as.matrix(predictors)

  S0 <- approx(model$info$blin_S$time, model$info$blin_S$S0, time)

  pred <- NULL
  for (row_number in 1:nrow(X)) {
    pred[row_number] <- S0$y[row_number]^as.numeric(exp(X[row_number, ]%*%model$beta + model$intercept))
  }

  predictions <-  cbind(pred, 1- pred)
  colnames(predictions) <- c("0", "1")
  return(predictions)
}

predict_linear_classification_numeric_survival <- function(model, predictors) {
  time <- predictors$time
  predictors$time <- NULL
  X <- as.matrix(predictors)

  as.numeric(X%*%model$beta + model$intercept)
}
