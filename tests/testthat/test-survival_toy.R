test_that("linear_classification() works", {
  set.seed(0)
  data <- simulate_dummy_surv_data()
  model <- linear_classification(event~time+., data)
  set.seed(1)
  new_data <- simulate_dummy_surv_data()
  pred <- predict(model, new_data)
  new_data$event
})

test_that("linear_classification() works with error_summary_CEN_surv()", {
  set.seed(0)
  data <- simulate_CEN_surv_data()
  model <- linear_classification(event~time+., data$data)
  set.seed(1)
  new_data <- simulate_CEN_surv_data()
  error_summary_CEN_surv(model, new_data)
})


test_that("linear_classification() works with parsnip", {
  library(parsnip)
  set.seed(0)
  data <- simulate_dummy_surv_data()
  model <- glasp_model() %>%
    set_mode("classification") %>%
    set_engine("glasp") %>%
    fit(event~., data)

  set.seed(1)
  new_data <- simulate_dummy_surv_data()
  pred <- predict(model, new_data, type = "prob")
  new_data$event
})

test_that("linear_classification() works with parsnip, tune and random search", {
  library(parsnip)
  library(tune)
  library(rsample)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_surv_data()
  model <- glasp_model(l1 = tune(), l2 = tune(), frob = tune(),
                       num_comp = tune()) %>%
    set_mode("classification") %>%
    set_engine("glasp")
  data_rs <- vfold_cv(data, v = 4)
  metric_vals <- metric_set(roc_auc)
  ctrl <- control_grid(verbose = T)

  hist <- tune_grid(model, event~.,
                        resamples = data_rs,
                        metrics = metric_vals,
                        control = ctrl,
                        grid = 10)
  show_best(hist, metric = "roc_auc")

  ctrl <- control_bayes(verbose = TRUE)
  hist <- tune_bayes(model, event~.,
                     resamples = data_rs,
                     metrics = metric_vals,
                     control = ctrl,
                     iter = 10)
  show_best(hist, metric = "roc_auc")
})

test_that("simulate_CEN_surv_data() works", {
  set.seed(0)
  data <- simulate_CEN_surv_data()
  model <- linear_classification(event~time+., data$data)
  set.seed(1)
  new_data <- simulate_CEN_surv_data()

  error_summary_CEN_surv(model, new_data)
})

test_that("Estimate baseline and survival function plot", {
  set.seed(1)

  N = 1000
  p = 10
  nzc = 5
  censProb = 0.3

  x = matrix(rnorm(N * p), N, p)
  beta = rnorm(nzc)
  fx = x[, seq(nzc)] %*% beta/3
  hx = exp(fx)
  ty = rexp(N, hx)
  tcens = rexp(N, hx/2)  # censoring indicator

  data <- data.frame(time = pmin(ty, tcens), event = factor((ty < tcens) + 0), x)

  model <- linear_classification(event~time + ., data, l1 = 0.03, l2 = 0.01, frob = 0.001, num_comp = 5)
  model$beta
  model$clusters

  N = 2
  x = matrix(rnorm(N * p), N, p)
  beta = rnorm(nzc)
  fx = x[, seq(nzc)] %*% beta/3
  hx = exp(fx)
  ty = rexp(N, hx)
  tcens = rexp(N, hx/2)  # censoring indicator
  new_data <- data.frame(time = pmin(ty, tcens), event = factor((ty < tcens) + 0), x)

  linear_predictor <- predict(model, new_data[1,], "numeric")
  time <- seq(0, 5, 0.1)
  S0 <- approx(model$info$blin_S$time, model$info$blin_S$S0, time)
  lambda <- hx[1]
  S_pred <- S0$y^exp(linear_predictor)
  S_true <- exp(-lambda*time)

  plot(S_pred)
  lines(S_true)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  data.frame(`GLASP estimation` = S_pred, `True function` = S_true, x = S0$x) %>% melt(id.vars = c("x")) %>%
  ggplot() + geom_line(aes(x = x, y = value, group = variable, linetype = variable, color = variable), size = 1) +
    #geom_point(aes(x = x, y = value, group = variable, color = variable)) +
    labs(color = "", linetype = "", y = "Survival function", x = "time") +
    theme_classic() +
    theme(legend.position = "top")

  })
