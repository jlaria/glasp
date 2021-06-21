#' Fit a `cox_regression`
#'
#' `cox_regression()` fits a model.
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
#' A `cox_regression` object.
#'
#' @examples
#' set.seed(0)
#' data <- simulate_dummy_surv_data()
#'
#' model <- cox_regression(event~time+., data, l1=0.05, l2=0.01, frob=0.01, num_comp=3)
#' model
#'

#' @export
cox_regression <- function(x, ...) {
  UseMethod("cox_regression")
}

#' @export
#' @rdname cox_regression
cox_regression.default <- function(x, ...) {
  stop("`cox_regression()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname cox_regression
cox_regression.data.frame <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  cox_regression_bridge(processed, ...)
}

# XY method - matrix

#' @export
#' @rdname cox_regression
cox_regression.matrix <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  cox_regression_bridge(processed, ...)
}

# Formula method

#' @export
#' @rdname cox_regression
cox_regression.formula <- function(formula, data, ...) {
  processed <- hardhat::mold(formula, data)
  cox_regression_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname cox_regression
cox_regression.recipe <- function(x, data, ...) {
  processed <- hardhat::mold(x, data)
  cox_regression_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

cox_regression_bridge <- function(processed, l1=0, l2=0, frob=0, num_comp=1, ...) {
  predictors <- processed$predictors
  #outcome <- processed$outcomes[[1]]
  outcome <- processed$outcomes


  hardhat::validate_predictors_are_numeric(predictors)
  hardhat::validate_outcomes_are_univariate(outcome)
  hardhat::validate_outcomes_are_binary(outcome)

  if(any(is.na(predictors)) | any(is.na(outcome))){
    stop("NA's are not supported")
  }

  hardhat::validate_column_names(predictors, "time")

  model <- cox_regression_impl(predictors, outcome, l1, l2, frob, num_comp)
  model$levels <- levels(outcome[[1]])


  new_cox_regression(
    model = model,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

cox_regression_impl <- function(predictors, outcome, l1, l2, frob, num_comp) {
  outcome[[1]] <- as.numeric(outcome[[1]]) - 1 # convert factor to 0-1 0:right censored, 1:event has been observed

  outcome <- cbind(predictors$time, outcome)
  predictors$time <- NULL

  dat <- list(X = as.matrix(predictors), y = as.matrix(outcome))
  obj <- new(glasp, dat, 2) # ph_cox
  obj$fit(c(l1, l2, frob, num_comp))

  beta <- as.numeric(obj$beta)
  names(beta) <- colnames(predictors)

  clusters = as.numeric(obj$clusters)
  names(clusters) = colnames(predictors)

  blin_S = compute_S0(dat$X, dat$y, obj$beta, obj$intercept)

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
        blin_S = blin_S
      ))
  )
}

compute_S0 <- function(X, y, beta, intercept){
  ord <- order(y[,1], decreasing = T)
  X <- X[ord, ]
  y <- y[ord, ]

  N <- nrow(X)

  timegrid <- y[, 1]
  eta <- intercept + X %*% beta
  S <- cumsum(exp(eta))

  h <- y[,2]/S

  timegrid <- rev(timegrid)
  #dif <- diff(c(0, timegrid))

  #H <- cumsum(rev(h)*dif)
  H <- cumsum(rev(h))

  S0 <- exp(-H)

  timegrid <- c(0, timegrid, Inf)
  S0 <- c(1, S0, 0)

  return(
    data.frame(time = timegrid,
               S0 = S0)
  )
}
