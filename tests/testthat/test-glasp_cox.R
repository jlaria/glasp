test_that("glasp_cox() fit and predict work", {
  library(glasp)
  library(parsnip)
  library(yardstick)

  set.seed(0)
  data <- simulate_dummy_surv_data(p=5, nzc=2)

  model <- glasp_cox(l1 = 0.005, l2 = 0.001, frob = 0.0001, num_comp = 3) %>%
    set_engine("glasp") %>% fit(event~time+., data)

  pred = predict(model$fit, data, type="class")
  pred$.truth = data$event
  acc = accuracy(pred, truth=.truth, estimate=.pred_class)

  expect_gt(acc$.estimate, 0.4)

  pred = predict(model$fit, data, type="prob")
  auc = roc_auc(pred, .pred_1, truth = data$event)
  expect_gt(auc$.estimate, 0.4)
})

test_that("glasp_cox() with tune works", {
  library(glasp)
  library(parsnip)
  library(tune)
  library(yardstick)
  library(rsample)

  set.seed(0)
  data <- simulate_dummy_surv_data(N=100, p=5, nzc=2)

  model <- glasp_classification(l1 = tune(),
                                l2 = tune(),
                                frob = tune(),
                                num_comp = tune()) %>% set_engine("glasp")

  data_rs <- vfold_cv(data, v = 4)

  hist <- tune_grid(model, event~time+.,
                    resamples = data_rs,
                    metrics = metric_set(roc_auc, accuracy),
                    grid = 5,
                    control = control_grid(verbose = F, save_pred = F)) #
  expect_gt(show_best(hist, 'roc_auc', 1)$mean, 0.4)
  expect_gt(show_best(hist, 'accuracy', 1)$mean, 0.4)
})
