test_that("cox_regression() fit works", {
  library(glasp)

  set.seed(0)
  data <- simulate_dummy_surv_data()
  model <- cox_regression(event~time+., data, l1=0, l2=0, frob=0, num_comp=1)

  expect_vector(model$model$beta)
})
