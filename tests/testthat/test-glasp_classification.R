test_that("glasp_classification() fit and predict work", {
  library(glasp)
  library(parsnip)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_logistic_data()

  model <- glasp_classification(l1 = 0.05, l2 = 0.01, frob = 0.001, num_comp = 3) %>%
    set_engine("glasp") %>% fit(y~., data)

  pred = predict(model$fit, data, type="class")
  pred$.truth = data$y
  acc = accuracy(pred, truth=.truth, estimate=.pred_class)

  expect_gt(acc$.estimate, 0.5)

  pred = predict(model$fit, data, type="prob")
  auc = roc_auc(pred, .pred_0, truth = data$y)
  expect_gt(auc$.estimate, 0.5)
})

test_that("glasp_classification() with tune works", {
  library(glasp)
  library(parsnip)
  library(tune)
  library(yardstick)
  library(rsample)

  set.seed(0)
  data <- simulate_dummy_logistic_data()

  model <- glasp_classification(l1 = tune(),
                                l2 = tune(),
                                frob = tune(),
                                num_comp = tune()) %>% set_engine("glasp")

  data_rs <- vfold_cv(data, v = 4)

  hist <- tune_grid(model, y~.,
                    resamples = data_rs,
                    metrics = metric_set(roc_auc, accuracy),
                    grid = 5,
                    control = control_grid(verbose = F, save_pred = F)) #

  expect_gt(show_best(hist, 'roc_auc', 1)$mean, 0.7)
  expect_gt(show_best(hist, 'accuracy', 1)$mean, 0.7)
})
