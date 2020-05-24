test_that("linear_regression() works", {
  library(parsnip)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_linear_data()
  model <- linear_regression(y~., data)
  set.seed(1)
  new_data <- simulate_dummy_linear_data()
  pred <- predict(model, new_data)
  rmse <- rmse_vec(new_data$y, pred)
  testthat::expect_lt(rmse, 1.55)

  model <- linear_regression(y~., data, l1 = 0.05)
  pred <- predict(model, new_data)
  rmse <- rmse_vec(new_data$y, pred)
  testthat::expect_lt(rmse, 1.27)
})

test_that("linear regression with parsnip works",{
  library(parsnip)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_linear_data()
  model <- glasp_model() %>%
    set_mode("regression") %>%
    set_engine("glasp") %>%
    fit(y~., data)
  set.seed(1)
  new_data <- simulate_dummy_linear_data()
  pred <- predict(model, new_data)
  rmse <- rmse_vec(new_data$y, pred$.pred)
  testthat::expect_lt(rmse, 1.55)

  model <- glasp_model(l1 = 0.05) %>%
    set_mode("regression") %>%
    set_engine("glasp") %>%
    fit(y~., data)
  pred <- predict(model, new_data)
  rmse <- rmse_vec(new_data$y, pred$.pred)
  testthat::expect_lt(rmse, 1.27)
})

test_that("NA's produce error", {
  set.seed(0)
  data <- simulate_dummy_linear_data()
  data[3,4]<- NA

  testthat::expect_error(model <- linear_regression(y~., data))
})
