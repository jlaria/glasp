test_that("linear_regression() fit and predict works", {
  library(glasp)

  set.seed(0)
  data <- simulate_dummy_linear_data()
  model <- linear_regression(y~., data, l1=0.01, l2=0.01, frob=0.001, num_comp=3)

  new_data <- simulate_dummy_linear_data()
  pred <- predict(model, new_data, type = "numeric")

  expect_vector(model$model$beta)
})
