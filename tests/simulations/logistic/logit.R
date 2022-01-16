library(glasp)
library(tune)
library(parsnip)
library(yardstick)
library(rsample)
library(doFuture)

all_cores <- parallel::detectCores(logical = FALSE) - 2
registerDoFuture()
cl <- parallel::makeCluster(all_cores)
plan(cluster, workers = cl)

results <- NULL

for (rho in c(0.8,0.5,0.2,0.1,0.0)) {
  for (p in c(100, 1000)) {
    for (run in 1:30) {
      sim <- simulate_CEN_logistic_data( N = 850, p = p, p_group = 50, true_p_group = 25, rho = rho)

      data_train <- sim$data[1:50,]
      data_test <- sim$data[-(1:50),]

      model <- glasp_classification(l1 = tune(), l2 = tune(), frob = tune(), num_comp = tune()) %>% set_engine('glasp')
      data_rs <- vfold_cv(data_train, v = 4)
      hist <- tune_grid(model, y~., data_rs, grid = 200, metrics = metric_set(accuracy), control = control_grid(verbose = F))

      best <- show_best(hist, n = 1)
      model <- glasp_classification(l1 = best$l1,
                                    l2 = best$l2,
                                    frob = best$frob,
                                    num_comp = best$num_comp) %>% set_engine('glasp') %>% fit(y~., data_train)


      acc <- accuracy_vec(data_test$y, predict(model, data_test, "class")[[1]])
      roc_auc <- roc_auc_vec(data_test$y, predict(model, data_test, "prob")[[1]])

      acc_beta <- mean((sim$true_beta!=0) == (abs(model$fit$model$beta) > 1e-4))
      num_iter <- length(model$fit$model$info$history)

      num_non_zeros <- sum(abs(model$fit$model$beta) > 1e-4)

      rand_index <- aricode::ARI(model$fit$model$clusters, sim$true_clusters_supervised)

      results <- rbind(results,
                       data.frame(
                         p = p,
                         rho = rho,
                         run = run,
                         acc = acc,
                         roc_auc = roc_auc,
                         correct_zeros = acc_beta,
                         num_non_zeros = num_non_zeros,
                         rand_index = rand_index,
                         num_iter = num_iter
                       ))
    }
  }
}
save.image('check.RData')
plan(sequential)

library(dplyr)

results %>% group_by(rho, p) %>% summarise(acc = mean(acc), roc_auc = mean(roc_auc), correct_zeros=mean(correct_zeros),
                                        rand_index = mean(rand_index), num_iter = mean(num_iter), num_non_zeros=mean(num_non_zeros))
