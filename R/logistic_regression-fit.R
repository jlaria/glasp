#' Fit a `logistic_regression`
#'
#' `logistic_regression()` fits a model.
#'
#' @param x Depending on the context:
#'
#'   * A __data frame__ of predictors.
#'   * A __matrix__ of predictors.
#'   * A __recipe__ specifying a set of preprocessing steps
#'     created from [recipes::recipe()].
#'
#' @param y When `x` is a __data frame__ or __matrix__, `y` is the outcome
#' specified as:
#'
#'   * A __data frame__ with 1 numeric column.
#'   * A __matrix__ with 1 numeric column.
#'   * A numeric __vector__.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both the predictors and the outcome.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#'
#' @return
#'
#' A `logistic_regression` object.
#'
#' @examples
#' set.seed(0)
#' data <- simulate_dummy_logistic_data()
#'
#' model <- logistic_regression(y~., data, l1=0.05, l2=0.01, frob=0.01, num_comp=3)
#' model
#'
#' new_data <- simulate_dummy_logistic_data()
#'
#' predict(model, new_data, type = "numeric")
#' predict(model, new_data, type = "prob")
#' predict(model, new_data, type = "class")

#' @export
logistic_regression <- function(x, ...) {
  UseMethod("logistic_regression")
}

#' @export
#' @rdname logistic_regression
logistic_regression.default <- function(x, ...) {
  stop("`logistic_regression()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname logistic_regression
logistic_regression.data.frame <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  logistic_regression_bridge(processed, ...)
}

# XY method - matrix

#' @export
#' @rdname logistic_regression
logistic_regression.matrix <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  logistic_regression_bridge(processed, ...)
}

# Formula method

#' @export
#' @rdname logistic_regression
logistic_regression.formula <- function(formula, data, ...) {
  processed <- hardhat::mold(formula, data)
  logistic_regression_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname logistic_regression
logistic_regression.recipe <- function(x, data, ...) {
  processed <- hardhat::mold(x, data)
  logistic_regression_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

logistic_regression_bridge <- function(processed, l1=0, l2=0, frob=0, num_comp=1, ...) {
  predictors <- processed$predictors
  #outcome <- processed$outcomes[[1]]
  outcome <- processed$outcomes


  hardhat::validate_predictors_are_numeric(predictors)
  hardhat::validate_outcomes_are_univariate(outcome)
  hardhat::validate_outcomes_are_binary(outcome)

  model <- logistic_regression_impl(predictors, outcome, l1, l2, frob, num_comp)
  model$levels <- levels(outcome[[1]])


  new_logistic_regression(
    model = model,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

logistic_regression_impl <- function(predictors, outcome, l1, l2, frob, num_comp) {

  outcome[[1]] <- as.numeric(outcome[[1]]) - 1
  dat <- list(X = as.matrix(predictors), y = as.matrix(outcome))
  obj <- new(glasp, dat, 1) # logistic
  obj$fit(c(l1, l2, frob, num_comp))

  beta <- as.numeric(obj$beta)
  names(beta) <- colnames(predictors)

  clusters = as.numeric(obj$clusters)
  names(clusters) = colnames(predictors)

  history = as.numeric(obj$history)

  return(
    list(
      beta = beta,
      intercept = obj$intercept,
      clusters = clusters,
      info = list(
        l1 = l1,
        l2 = l2,
        frob = frob,
        num_comp = num_comp,
        history = history
      ))
    )
}

