#' Simulate toy linear data for benchmark
#'
#' @param N Number of observations
#' @param p Number of features
#' @param nzc True number of features in the model
#'
#' @return An expanded data set
#' @export
#'
#' @importFrom "stats" "rnorm"
simulate_dummy_linear_data <- function(N = 100, p = 10, nzc = 5){
  x = matrix(rnorm(N * p), N, p)
  beta = rnorm(nzc)
  fx = x[, seq(nzc)] %*% beta/3
  y = fx + rnorm(N)

  data <- data.frame(y = y, x)
}

#' Simulate toy logistic data for benchmark
#'
#' @param N Number of observations
#' @param p Number of features
#' @param nzc True number of features in the model
#'
#' @return An expanded data set
#' @export
#'
#' @importFrom "stats" "rnorm" "rbinom"
simulate_dummy_logistic_data <- function(N = 100, p = 10, nzc = 5){
  x = matrix(rnorm(N * p), N, p)
  beta = rnorm(nzc)
  fx = x[, seq(nzc)] %*% beta
  eta = fx + rnorm(N)
  p = (1 + exp(-eta))^-1
  y = rbinom(N, 1, p)
  data <- data.frame(y = factor(y, levels = c('0', '1')), x)
}


rmse_beta <- function(true, pred, Sigma){
  #sqrt(mean((true - pred)^2))
  t(true - pred)%*%Sigma%*%(true-pred)
}


#' Simulate toy survival data for benchmark
#'
#' @param N Number of observations (event or censored time)
#' @param p Number of features
#' @param nzc True number of features in the model
#' @param censProb Probability of right censoring
#'
#' @return An expanded data set
#' @export
#'
#' @importFrom "stats" "rexp"
#' @importFrom "stats" "rnorm"
simulate_dummy_surv_data <- function(N = 100, p = 10, nzc = 5, censProb = 0.3){
  x = matrix(rnorm(N * p), N, p)
  beta = rnorm(nzc)
  fx = x[, seq(nzc)] %*% beta/3
  hx = exp(fx)
  ty = rexp(N, hx)
  tcens = rexp(N, hx/2)  # censoring indicator

  data <- data.frame(time = pmin(ty, tcens), event = factor((ty < tcens) + 0), x)
}

#' Simulate survival data with groups
#'
#' @param N Number of observations (event or censored time)
#' @param p Number of features
#' @param p_group Number of features in each group
#' @param true_p_group Number of important features in each group
#' @param rho Correlation inside groups
#'
#' @return A data set0   2
#' @export
#'
#' @importFrom "stats" "rexp"
#' @importFrom "MASS" "mvrnorm"
#' @importFrom "stats" "rnorm"
#' @importFrom "stats" "runif"
#'
simulate_CEN_surv_data <- function(N = 100, p = 30, p_group = 10, true_p_group = 3, rho = 0.5){
  sigma_star <- matrix(rho, nrow = p_group, ncol = p_group)
  diag(sigma_star) <- 1

  X_1 <- mvrnorm(N, mu = rep(0, p_group), Sigma = sigma_star)
  X_2 <- mvrnorm(N, mu = rep(0, p_group), Sigma = sigma_star)
  X_3 <- matrix(rnorm((p-2*p_group)*N), nrow = N, ncol = p - 2*p_group)

  x <- cbind(X_1, X_2, X_3)
  beta <- c(
    runif(true_p_group, 0.9, 1.1),
    rep(0, p_group - true_p_group),
    -runif(true_p_group, 0.9, 1.1),
    rep(0, p - p_group - true_p_group)
  )

  fx = x %*% beta
  hx = exp(fx)
  ty = rexp(N, hx)
  tcens = rexp(N, hx/2)  # censoring indicator

  data <- data.frame(time = pmin(ty, tcens), event = factor((ty < tcens) + 0), x)

  Sigma <- diag(nrow = p, ncol = p)
  Sigma[1:p_group, 1:p_group] <- sigma_star
  Sigma[(p_group + 1):(2*p_group), (p_group + 1):(2*p_group)] <- sigma_star

  list(
    data = data,
    true_beta = beta,
    true_clusters = rep(1:3, times = c(p_group, p_group, p - 2*p_group)),
    true_clusters_supervised = c(rep(1, true_p_group), rep(2, p_group - true_p_group),
                                 rep(3, true_p_group), rep(2, p_group - true_p_group),
                                 rep(2, p - 2*p_group)),
    lambdas = hx,
    Sigma = Sigma
  )
}

#' Simulate classification data with groups
#'
#' @param N Number of observations (event or censored time)
#' @param p Number of features
#' @param p_group Number of features in each group
#' @param true_p_group Number of important features in each group
#' @param rho Correlation inside groups
#'
#' @export
#'
#' @importFrom "MASS" "mvrnorm"
#'
simulate_CEN_logistic_data <- function(N = 100, p = 30, p_group = 10, true_p_group = 3, rho = 0.5){
  sigma_star <- matrix(rho, nrow = p_group, ncol = p_group)
  diag(sigma_star) <- 1

  X_1 <- mvrnorm(N, mu = rep(0, p_group), Sigma = sigma_star)
  X_2 <- mvrnorm(N, mu = rep(0, p_group), Sigma = sigma_star)
  X_3 <- matrix(rnorm((p-2*p_group)*N), nrow = N, ncol = p - 2*p_group)

  x <- cbind(X_1, X_2, X_3)
  beta <- c(
    runif(true_p_group, 0.9, 1.1),
    rep(0, p_group - true_p_group),
    -runif(true_p_group, 0.9, 1.1),
    rep(0, p - p_group - true_p_group)
  )

  eta = x %*% beta
  prob = (1+exp(-eta))^-1
  y = factor((prob > 0.5) + 0)
  data <- data.frame(y = y, x)

  Sigma <- diag(nrow = p, ncol = p)
  Sigma[1:p_group, 1:p_group] <- sigma_star
  Sigma[(p_group + 1):(2*p_group), (p_group + 1):(2*p_group)] <- sigma_star

  list(
    data = data,
    true_beta = beta,
    true_clusters = rep(1:3, times = c(p_group, p_group, p - 2*p_group)),
    true_clusters_supervised = c(rep(1, true_p_group), rep(2, p_group - true_p_group),
                                 rep(3, true_p_group), rep(2, p_group - true_p_group),
                                 rep(2, p - 2*p_group)),
    Sigma = Sigma
  )
}

