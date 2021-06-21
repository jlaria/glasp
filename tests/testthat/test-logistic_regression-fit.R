test_that("logistic_regression() fit works", {
  library(glasp)

  set.seed(0)
  data <- simulate_dummy_logistic_data()
  model <- logistic_regression(y~., data, l1=0, l2=0, frob=0, num_comp=1)

  expect_vector(model$model$beta)
})

