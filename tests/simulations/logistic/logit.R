library(glasp)
library(tune)
library(parsnip)
library(yardstick)
library(rsample)
library(doFuture)
library(glue)

all_cores <- 6
registerDoFuture()
cl <- parallel::makeCluster(all_cores)
plan(cluster, workers = cl)

p = 1000
for (rho in c(0.8,0.5,0.2,0.1,0.0)) {
    for (run in 1:30) {
      sim <- simulate_CEN_logistic_data( N = 850, p = p, p_group = 50, true_p_group = 25, rho = rho)

      data_train <- sim$data[1:200,]
      data_test <- sim$data[-(1:200),]

      model <- glasp_classification(l1 = tune(), l2 = tune(), frob = tune(), num_comp = tune()) %>% set_engine('glasp')
      data_rs <- vfold_cv(data_train, v = 2)
      hist <- tune_grid(model, y~., data_rs, grid = 100, metrics = metric_set(accuracy), control = control_grid(verbose = F))

      best <- show_best(hist, n = 1)
      model <- glasp_classification(l1 = best$l1,
                                    l2 = best$l2,
                                    frob = best$frob,
                                    num_comp = best$num_comp) %>% set_engine('glasp') %>% fit(y~., data_train)


      acc <- accuracy_vec(data_test$y, predict(model, data_test, "class")[[1]])
      roc_auc <- roc_auc_vec(data_test$y, predict(model, data_test, "prob")[[1]])

      acc_beta <- mean((sim$true_beta!=0) == (abs(model$fit$model$beta) > 1e-2))
      num_iter <- length(model$fit$model$info$history)

      num_non_zeros <- sum(abs(model$fit$model$beta) > 1e-2)

      rand_index_supervised <- aricode::ARI(model$fit$model$clusters, sim$true_clusters_supervised)
      rand_index_unsupervised <- aricode::ARI(model$fit$model$clusters, sim$true_clusters)

      write.csv(
      data.frame(p = p,
                 rho = rho,
                 run = run,
                 acc = acc,
                 roc_auc = roc_auc,
                 correct_zeros = acc_beta,
                 num_non_zeros = num_non_zeros,
                 rand_index_supervised = rand_index_supervised,
                 rand_index_unsupervised = rand_index_unsupervised,
                 num_iter = num_iter
               ),
      file=glue('/workspace/code/glasp/tests/logistic/logs/{rho}-{p}-{run}.csv'))
      print(glue('logs/{rho}-{p}-{run}.csv'))
    }
  }
print('Finished!')

plan(sequential)

library(dplyr)

rm(list=ls())

files <- dir('/workspace/code/glasp/tests/logistic/logs/')
results <- NULL

for( f in files ){
  in_ <- read.csv(paste0('/workspace/code/glasp/tests/logistic/logs/',f))
  results <- rbind(results, in_)
}



results %>% group_by(rho, p) %>% summarise(acc = mean(acc), roc_auc = mean(roc_auc), correct_zeros=mean(correct_zeros),
                                        rand_index_supervised = mean(rand_index_supervised),
                                        rand_index_unsupervised = mean(rand_index_unsupervised),
                                        num_iter = mean(num_iter), num_non_zeros=mean(num_non_zeros))
# %>%  write.csv(file='/workspace/code/glasp/tests/logistic/results.csv')
results %>% group_by(rho, p) %>% summarise(acc = sd(acc)/sqrt(30), roc_auc = sd(roc_auc)/sqrt(30),
                                           correct_zeros=sd(correct_zeros)/sqrt(30),
                                           rand_index_supervised = sd(rand_index_supervised)/sqrt(30),
                                           rand_index_unsupervised = sd(rand_index_unsupervised)/sqrt(30),
                                           num_iter = sd(num_iter)/sqrt(30),
                                           num_non_zeros=sd(num_non_zeros)/sqrt(30))
