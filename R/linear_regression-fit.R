#' Fit a `linear_regression`
#'
#' `linear_regression()` fits a model.
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
#' A `linear_regression` object.
#'
#' @examples
#' set.seed(0)
#' data <- simulate_dummy_linear_data()
#'
#' model <- linear_regression(y~., data, l1=0.05, l2=0.01, frob=0.01, num_comp=3)
#' model
#'
#' new_data <- simulate_dummy_linear_data()
#'
#' predict(model, new_data, type = "numeric")

#' @export
linear_regression <- function(x, ...) {
  UseMethod("linear_regression")
}

#' @export
#' @rdname linear_regression
linear_regression.default <- function(x, ...) {
  stop("`linear_regression()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname linear_regression
linear_regression.data.frame <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  linear_regression_bridge(processed, ...)
}

# XY method - matrix

#' @export
#' @rdname linear_regression
linear_regression.matrix <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  linear_regression_bridge(processed, ...)
}

# Formula method

#' @export
#' @rdname linear_regression
linear_regression.formula <- function(formula, data, ...) {
  processed <- hardhat::mold(formula, data)
  linear_regression_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname linear_regression
linear_regression.recipe <- function(x, data, ...) {
  processed <- hardhat::mold(x, data)
  linear_regression_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

linear_regression_bridge <- function(processed, l1=0, l2=0, frob=0, num_comp=1, ...) {
  predictors <- processed$predictors
  #outcome <- processed$outcomes[[1]]
  outcome <- processed$outcomes


  hardhat::validate_predictors_are_numeric(predictors)
  hardhat::validate_outcomes_are_univariate(outcome)
  hardhat::validate_outcomes_are_numeric(outcome)

  model <- linear_regression_impl(predictors, outcome, l1, l2, frob, num_comp)

  new_linear_regression(
    model = model,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

linear_regression_impl <- function(predictors, outcome, l1, l2, frob, num_comp) {

  dat <- list(X = as.matrix(predictors), y = as.matrix(outcome))
  obj <- new(glasp, dat, 0) # linear
  obj$fit(c(l1, l2, frob, num_comp))

  beta <- as.numeric(obj$beta)
  names(beta) <- colnames(predictors)

  clusters = as.numeric(obj$clusters)
  names(clusters) = colnames(predictors)


  return(
    list(
      beta = beta,
      intercept = obj$intercept,
      clusters = clusters,
      info = list(
        l1 = l1,
        l2 = l2,
        frob = frob,
        num_comp = num_comp
      ))
  )
}

