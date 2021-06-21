test_that("cox_regression() fit and predict work", {
  library(glasp)

  set.seed(0)
  data <- simulate_dummy_surv_data(N=1000)
  model <- cox_regression(event~time+., data, l1=0.01, l2=0.01, frob=0.001, num_comp=3)

  new_data <- simulate_dummy_surv_data()
  pred <- predict(model, new_data, type = "numeric")
  pred <- predict(model, new_data, type = "prob")
  pred <- predict(model, new_data, type = "class")

  expect_vector(model$model$beta)
})
