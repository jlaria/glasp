test_that("logistic_regression() fit and predict works", {
  library(glasp)

  set.seed(0)
  data <- simulate_dummy_logistic_data()
  model <- logistic_regression(y~., data, l1=0.1, l2=0.01, frob=0.001, num_comp=3)
  #model$model$beta
  #model$model$info$history

  new_data <- simulate_dummy_logistic_data()
  pred <- predict(model, new_data, type = "numeric")
  pred <- predict(model, new_data, type = "prob")
  pred <- predict(model, new_data, type = "class")

  expect_vector(model$model$beta)
})
