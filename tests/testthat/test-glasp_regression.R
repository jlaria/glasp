test_that("glasp_regression() fit and predict work", {
  library(glasp)
  library(parsnip)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_linear_data()

  model <- glasp_regression(l1 = 0.05, l2 = 0.01, frob = 0.001, num_comp = 3) %>%
    set_engine("glasp") %>% fit(y~., data)

  pred = predict(model$fit, data, type="numeric")
  rmse = rmse_vec(data$y, pred)
  expect_lt(rmse, 1.5)
})

test_that("glasp_regression() with tune works", {
  library(glasp)
  library(parsnip)
  library(tune)
  library(yardstick)
  library(rsample)

  set.seed(0)
  data <- simulate_dummy_linear_data()

  model <- glasp_regression(    l1 = tune(),
                                l2 = tune(),
                                frob = tune(),
                                num_comp = tune()) %>% set_engine("glasp")

  data_rs <- vfold_cv(data, v = 4)

  hist <- tune_grid(model, y~.,
                    resamples = data_rs,
                    metrics = metric_set(yardstick::rmse),
                    grid = 5,
                    control = control_grid(verbose = F, save_pred = F)) #

  best <- select_best(hist, 'rmse', 1)
  final_model <- glasp_regression(l1 = best$l1, l2 = best$l2, frob = best$frob, num_comp = best$num_comp) %>%
    set_engine("glasp") %>% fit(y~., data)

  expect_lt(rmse_vec(data$y, predict(final_model$fit, data)), 2.0)
})
