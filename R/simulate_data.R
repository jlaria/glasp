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
  data <- data.frame(y = factor(y), x)
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


#' Evaluate survival simulations
#'
#'
#' @param model A `linear_classification` model
#' @param new_data Data generated from `simulate_CEN_surv_data`
#'
#' @importFrom "fda.usc" "fdata"
#' @importFrom "fda.usc" "int.simpson"
#' @importFrom "aricode" "ARI"
#' @export
error_summary_CEN_surv <- function(model, new_data){
  # 1 - beta error
  beta_rmse <- rmse_beta(new_data$true_beta, model$beta, new_data$Sigma)

  # beta ccr
  beta_numnonzeros <- sum(model$beta!=0)
  beta_correctzeros <- mean((model$beta!=0) == (new_data$true_beta!=0))
  beta_tpr <- mean((model$beta!=0)[new_data$true_beta!=0])
  beta_tnr <- mean((model$beta==0)[new_data$true_beta==0])
  # 2 - cluster info
  non_zero_coef <- new_data$true_beta != 0

#  RI <- adjustedRand(new_data$true_clusters, model$clusters, randMethod="Rand")
  RI <- ARI(new_data$true_clusters, model$clusters)
#  RI_sup <- adjustedRand(new_data$true_clusters_supervised[non_zero_coef], model$clusters[non_zero_coef], randMethod = "Rand")
  RI_sup <- ARI(new_data$true_clusters_supervised[non_zero_coef], model$clusters[non_zero_coef])
  vi.dist <- vi_dist(new_data$true_clusters, model$clusters)
  vi.dist_sup <- vi_dist(new_data$true_clusters_supervised[non_zero_coef], model$clusters[non_zero_coef])

  # 3 -  S error
  eta <- predict.linear_classification(model, new_data$data, "numeric")
  time <- sort(new_data$data$time)
  S0 <- approx(model$info$blin_S$time, model$info$blin_S$S0, time)

  abs_S_error <- NULL
  for (i in 1:length(eta)) {
    linear_predictor <- eta[i]
    lambda <- new_data$lambdas[i]

    S_pred <- S0$y^exp(linear_predictor)
    S_true <- exp(-lambda*time)

    f <- fdata(abs(S_pred - S_true))
    abs_S_error[i] <- int.simpson(f)
  }

  mean_abs_S_error <- mean(abs_S_error)

  data.frame(
    beta_rmse = beta_rmse,
    beta_correctzeros = beta_correctzeros,
    beta_numnonzeros = beta_numnonzeros,
    vi.dist = vi.dist,
    vi.dist_sup = vi.dist_sup,
    mean_abs_S_error = mean_abs_S_error,
    beta_tpr = beta_tpr,
    beta_tnr = beta_tnr,
    RI = RI,
    RI_sup = RI_sup
  )
}
