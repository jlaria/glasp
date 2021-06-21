#' Predict from a `logistic_regression`
#'
#' @param object A `logistic_regression` object.
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
#' @examples
#' set.seed(0)
#' data <- simulate_dummy_logistic_data()
#'
#' model <- logistic_regression(y~., data, l1=0.05, l2=0.01, frob=0.01, num_comp=3)
#'
#' new_data <- simulate_dummy_logistic_data()
#'
#' predict(model, new_data, type = "numeric")
#' predict(model, new_data, type = "prob")
#' predict(model, new_data, type = "class")
#'
#' @export
predict.logistic_regression <- function(object, new_data, type = "numeric", ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_logistic_regression_predict_types())
  predict_logistic_regression_bridge(type, object, forged$predictors)
}

valid_logistic_regression_predict_types <- function() {
  c("numeric", "prob", "class")
}

# ------------------------------------------------------------------------------
# Bridge

predict_logistic_regression_bridge <- function(type, model, predictors) {
  #predictors <- as.matrix(predictors)

  predict_function <- get_logistic_regression_predict_function(type)
  predictions <- predict_function(model, predictors)

  hardhat::validate_prediction_size(predictions, predictors)

  predictions
}

get_logistic_regression_predict_function <- function(type) {
  switch(
    type,
    numeric = predict_logistic_regression_numeric,
    prob = predict_logistic_regression_prob,
    class = predict_logistic_regression_class,
  )
}

# ------------------------------------------------------------------------------
# Implementation

predict_logistic_regression_numeric <- function(model, predictors) {
  X <- as.matrix(predictors)
  eta <- as.numeric(X%*%model$model$beta + model$model$intercept)

  #hardhat::spruce_numeric(as.numeric((1+exp(-eta))^-1))
  (as.numeric((1+exp(-eta))^-1) > 0.5) + 0
}

predict_logistic_regression_prob <- function(model, predictors) {
  X <- as.matrix(predictors)
  eta <- as.numeric(X%*%model$model$beta + model$model$intercept)

  pred <- as.numeric((1+exp(-eta))^-1)
  predictions <-  cbind(1 - pred, pred)

  hardhat::spruce_prob(pred_levels = model$model$levels, prob_matrix = predictions)
}

predict_logistic_regression_class <- function(model, predictors) {
  X <- as.matrix(predictors)
  eta <- as.numeric(X%*%model$model$beta + model$model$intercept)

  pred <- (as.numeric((1+exp(-eta))^-1) > 0.5) + 1
  hardhat::spruce_class(pred_class = factor(model$model$levels[pred]))
}
