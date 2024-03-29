% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glasp_classification.R
\name{glasp_classification}
\alias{glasp_classification}
\title{parsnip interface for GLASP - Classification Models}
\usage{
glasp_classification(
  mode = "classification",
  l1 = NULL,
  l2 = NULL,
  frob = NULL,
  num_comp = NULL
)
}
\arguments{
\item{mode}{A single character string for the type of model.}

\item{l1}{Regularization in norm l1 (lasso)}

\item{l2}{Regularization in norm l2 (ridge)}

\item{frob}{Regularization in Frobenius norm. It is unique to glasp and controls the importance
of the clustering on the variables. If frob=0, we want to fit a model equivalent to elastic-net,
without clustering. If frob > l1+l2 means that we are more interested on the feature clustering
than the variable selection. If frob > 1 we are more interested on finding clusters than on
getting a linear model with good predictive power.}

\item{num_comp}{Maximum number of clusters to search.}
}
\value{
A \code{glasp_classification} parsnip model
}
\description{
parsnip interface for GLASP - Classification Models
}
\examples{

library(parsnip)
library(yardstick)

set.seed(0)
data <- simulate_dummy_logistic_data()

model <- glasp_classification(l1 = 0.05, l2 = 0.01, frob = 0.001, num_comp = 3) \%>\%
  set_engine("glasp") \%>\%
  fit(y~., data)

pred = predict(model$fit, data, type="class")
pred$.truth = data$y
acc = accuracy(pred, truth=.truth, estimate=.pred_class)

pred = predict(model$fit, data, type="prob")
auc = roc_auc(pred, .pred_0, truth = data$y)

#--- with tune ----
\dontrun{
library(parsnip)
library(tune)
library(yardstick)
library(rsample)

set.seed(0)
data <- simulate_dummy_logistic_data()

model <- glasp_classification(l1 = tune(), l2 = tune(), frob = tune(), num_comp = tune()) \%>\%
set_engine("glasp")

data_rs <- vfold_cv(data, v = 4)
hist <- tune_grid(model, y~.,
        resamples = data_rs,
        metrics = metric_set(roc_auc, accuracy),
        grid = 100,
        control = control_grid(verbose = T, save_pred = T))
show_best(hist, 'roc_auc')
show_best(hist, 'accuracy')
}

}
