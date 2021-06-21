#' Predict from a `cox_regression`
#'
#' @param object A `cox_regression` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"numeric"` for numeric predictions.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @export
predict.cox_regression <- function(object, new_data, type = "numeric", ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_cox_regression_predict_types())
  predict_cox_regression_bridge(type, object, forged$predictors)
}

valid_cox_regression_predict_types <- function() {
  c("numeric", "prob", "class")
}

# ------------------------------------------------------------------------------
# Bridge

predict_cox_regression_bridge <- function(type, model, predictors) {
  #predictors <- as.matrix(predictors)

  predict_function <- get_cox_regression_predict_function(type)
  predictions <- predict_function(model, predictors)

  hardhat::validate_prediction_size(predictions, predictors)

  predictions
}

get_cox_regression_predict_function <- function(type) {
  switch(
    type,
    numeric = predict_cox_regression_numeric,
    prob = predict_cox_regression_prob,
    class = predict_cox_regression_class,
  )
}

# ------------------------------------------------------------------------------
# Implementation

#' @importFrom "stats" "approx"
get_pred <- function(model, predictors){
  # 1) We need the baseline survival function model$blin_S
  # 2) predictors have to be in the expanded format

  time <- predictors$time

  predictors$time <- NULL

  X <- as.matrix(predictors)

  S0 <- approx(model$model$info$blin_S$time, model$model$info$blin_S$S0, time)

  pred <- NULL
  for (row_number in 1:nrow(X)) {
    pred[row_number] <- S0$y[row_number]^as.numeric(exp(X[row_number, ]%*%model$model$beta + model$model$intercept))
  }
  return (pred)
}


predict_cox_regression_numeric <- function(model, predictors) {
  pred <- get_pred(model, predictors)
  as.numeric((pred > 0.5) + 0)
}

predict_cox_regression_prob <- function(model, predictors) {

  pred <- get_pred(model, predictors)
  predictions <-  cbind(1 - pred, pred)

  hardhat::spruce_prob(pred_levels = model$model$levels, prob_matrix = predictions)
}

predict_cox_regression_class <- function(model, predictors) {
  pred <- as.numeric(get_pred(model, predictors) > 0.5) + 1
  hardhat::spruce_class(pred_class = factor(model$model$levels[pred]))
}
