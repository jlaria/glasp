
<!-- badges: start -->

[![R-CMD-check](https://github.com/jlaria/glasp/workflows/R-CMD-check/badge.svg)](https://github.com/jlaria/glasp/actions)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Group Linear Algorithm with Sparse Principal decomposition <img src="man/figures/glasp.png" align="right" width="150" />

> A variable selection and clustering method for generalized linear
> models

## Overview

R package `glasp` is an extension of the Sparse Group Lasso that
computes groups automatically. The internal supervised variable
clustering algorithm is also an original contribution and integrates
naturally within the Group Lasso penalty. Moreover, this implementation
provides the flexibility to change the risk function and address any
regression problem.

## Install

To install `glasp` in R use

    devtools::install_github("jlaria/glasp", dependencies = TRUE)

## Example

### Linear regression

In this example, we show how to integrate `glasp` with `parsnip` and
other tidy libraries.

We first load the required libraries.

``` r
library(glasp)
library(parsnip)
library(yardstick)
#> For binary classification, the first factor level is assumed to be the event.
#> Use the argument `event_level = "second"` to alter this as needed.
```

Next, we simulate some linear data with the `simulate_dummy_linear_data`
function.

``` r
set.seed(0)
data <- simulate_dummy_linear_data()
```

A `glasp` model can be computed using different approaches. This is the
`parsnip` approach.

``` r
model <- glasp_regression(l1=0.01, l2=0.001, frob=0.0003) %>%
    set_engine("glasp") %>%
    fit(y~., data)
```

The coefficients can be accessed through the `parsnip` model object
`model$fit$beta`.

``` r
print(model)
#> parsnip model object
#> 
#> Fit time:  25ms 
#> <linear_regression> 
#> $model
#> $model$beta
#>           X1           X2           X3           X4           X5           X6 
#>  0.004520053  0.402060668 -0.032605834 -0.369015410 -0.689548268  0.000000000 
#>           X7           X8           X9          X10 
#>  0.000000000  0.023018621  0.000000000  0.000000000 
#> 
#> $model$intercept
#> [1] -0.1668452
#> 
#> $model$clusters
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   1   1   1   1   0   1   1   1   1   1 
#> 
#> $model$info
#> $model$info$l1
#> [1] 0.01
#> 
#> $model$info$l2
#> [1] 0.001
#> 
#> $model$info$frob
#> [1] 3e-04
#> 
#> $model$info$num_comp
#> [1] 1
```

Now we generate some out-of-sample data to check how the model predicts.

``` r
set.seed(1)
new_data <- simulate_dummy_linear_data()

pred <- predict(model, new_data)
rmse <- rmse_vec(new_data$y, pred$.pred)
```

We obtain a 1.5180255 root mean square error.

### Hyper-parameter selection

Hyper-parameter search is quite easy using the tools `glasp` integrates
with. To show this, we will simulate some survival data.

``` r
set.seed(0)
data <- simulate_dummy_surv_data()
head(data)
#>         time event         X1         X2         X3         X4         X5
#> 1 0.37444362     1  1.2629543  0.7818592 -1.0457177 -0.1244350 -0.3900096
#> 2 0.70084082     1 -0.3262334 -0.7767766 -0.8962113  1.4667446 -1.8192222
#> 3 0.66139233     1  1.3297993 -0.6159899  1.2693872  0.6739287  0.6591807
#> 4 2.14750885     0  1.2724293  0.0465803  0.5938409  1.9564253  0.4596217
#> 5 4.31371385     0  0.4146414 -1.1303858  0.7756343 -0.2690410  1.6166263
#> 6 0.05265771     0 -1.5399500  0.5767188  1.5573704 -1.2445515 -1.8561905
#>           X6         X7         X8         X9         X10
#> 1 -0.6263682 -1.6878010  0.5812010 -1.4240317 -0.59188422
#> 2  0.4813353  0.6476460  0.3792931 -1.6693444 -0.37099306
#> 3  1.6952711  0.4487942 -0.3107409  1.3792361  0.08792426
#> 4 -1.7612263  1.0263022  0.8863900 -0.9196746 -0.03472634
#> 5  0.1980130  1.0749782 -1.6418647 -0.5044900  1.80637427
#> 6  0.3973491  0.4583096 -0.9885637 -1.1347318 -0.34023607
```

We create the glasp model, but this time we do not fit it. Instead, we
will call the `tune` function.

``` r
library(tune)
#> Registered S3 method overwritten by 'tune':
#>   method                   from   
#>   required_pkgs.model_spec parsnip

model <- glasp_cox(l1 = tune(), 
                     l2 = tune(),
                     frob = tune(), 
                     num_comp = tune()) %>%
         set_engine("glasp")
```

We specify k-fold cross validation and Bayesian optimization to search
for hyper-parameters.

``` r
library(rsample)

data_rs <- vfold_cv(data, v = 4)

hist <- tune_bayes(model, event~time+., # <- Notice the syntax with time in the right hand side 
                   resamples = data_rs,
                   metrics = metric_set(roc_auc), # yardstick's roc_auc
                   iter =10, # <- 10 iterations... change to 1000
                   control = control_bayes(verbose = TRUE)) # 
#> 
#> >  Generating a set of 5 initial parameter results
#> ✓ Initialization complete
#> 
#> Optimizing roc_auc using the expected improvement
#> 
#> ── Iteration 1 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.558 (@iter 0)
#> i Gaussian process model
#> ! The Gaussian process model is being fit using 4 features but only has 5
#>   data points to do so. This may cause errors or a poor model fit.
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=6.17e-06, l2=1.88e-06, frob=0.00377, num_comp=9
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.5606 (+/-0.0787)
#> 
#> ── Iteration 2 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=2.55, l2=3.57, frob=0.0911, num_comp=17
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 3 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.00564, l2=0.000567, frob=9.93, num_comp=18
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5587 (+/-0.0599)
#> 
#> ── Iteration 4 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.000125, l2=2.25, frob=0.0229, num_comp=19
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 5 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.0128, l2=3.38e-06, frob=0.0265, num_comp=7
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.548 (+/-0.0737)
#> 
#> ── Iteration 6 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.0497, l2=0.000154, frob=4.47, num_comp=9
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5603 (+/-0.0602)
#> 
#> ── Iteration 7 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=8.6e-06, l2=3.57, frob=1.12e-06, num_comp=10
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 8 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.018, l2=4.1e-05, frob=3.7, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5587 (+/-0.0606)
#> 
#> ── Iteration 9 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=3.55, l2=2.04e-06, frob=6.97e-06, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 10 ────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5606 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.000668, l2=5.35e-06, frob=4.39e-06, num_comp=10
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5382 (+/-0.0825)

show_best(hist, metric = "roc_auc")
#> # A tibble: 5 × 11
#>           l1         l2    frob num_comp .metric .estimator  mean     n std_err
#>        <dbl>      <dbl>   <dbl>    <int> <chr>   <chr>      <dbl> <int>   <dbl>
#> 1 0.00000617 0.00000188 0.00377        9 roc_auc binary     0.561     4  0.0787
#> 2 0.0497     0.000154   4.47           9 roc_auc binary     0.560     4  0.0602
#> 3 0.00564    0.000567   9.93          18 roc_auc binary     0.559     4  0.0599
#> 4 0.0180     0.0000410  3.70           4 roc_auc binary     0.559     4  0.0606
#> 5 0.0000248  0.00000669 0.122         17 roc_auc binary     0.558     4  0.0683
#> # … with 2 more variables: .config <chr>, .iter <int>
```

### Logistic regression

``` r
library(glasp)
library(yardstick)

set.seed(0)
data <- simulate_dummy_logistic_data()
model <- logistic_regression(y~., data, l1=0.01, l2=0.001, frob=0.001, ncomp=2)
print(model)
#> <logistic_regression> 
#> $model
#> $model$beta
#>          X1          X2          X3          X4          X5          X6 
#> -0.11957419  0.51303973 -0.13348002 -0.44564217 -1.14015520  0.10695661 
#>          X7          X8          X9         X10 
#>  0.07211584  0.02068296 -0.04030295  0.00000000 
#> 
#> $model$intercept
#> [1] -0.01230832
#> 
#> $model$clusters
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   1   1   1   1   0   1   1   1   1   1 
#> 
#> $model$info
#> $model$info$l1
#> [1] 0.01
#> 
#> $model$info$l2
#> [1] 0.001
#> 
#> $model$info$frob
#> [1] 0.001
#> 
#> $model$info$num_comp
#> [1] 1
#> 
#> $model$info$history
#> [1] 2.237973e-01 5.495960e-02 8.738682e-03 1.553462e-03 2.813430e-04
#> [6] 4.522578e-05
#> 
#> 
#> $model$levels
#> [1] "0" "1"

pred = predict(model, data)
```

``` r
library(rsample)

set.seed(0)
data <- simulate_dummy_logistic_data()

model <- glasp_classification(l1 = tune(), 
                     l2 = tune(),
                     frob = tune(), 
                     num_comp = tune()) %>%
         set_engine("glasp")

data_rs <- vfold_cv(data, v = 4)

hist <- tune_grid(model, y~., 
                   resamples = data_rs,
                   metrics = metric_set(roc_auc), 
                   grid =10, 
                   control = control_grid(verbose = FALSE)) # 
show_best(hist, metric = "roc_auc")
#> # A tibble: 5 × 10
#>           l1         l2     frob num_comp .metric .estimator  mean     n std_err
#>        <dbl>      <dbl>    <dbl>    <int> <chr>   <chr>      <dbl> <int>   <dbl>
#> 1 0.000520   0.00841     7.01e-4       17 roc_auc binary     0.865     4  0.0214
#> 2 0.0161     0.000998    1.84e-4       19 roc_auc binary     0.862     4  0.0237
#> 3 0.00000192 0.00000327  1.69e-6       10 roc_auc binary     0.860     4  0.0255
#> 4 0.0000258  0.0312      7.34e-6       12 roc_auc binary     0.859     4  0.0237
#> 5 0.00795    0.00000599  1.88e-2       15 roc_auc binary     0.833     4  0.0223
#> # … with 1 more variable: .config <chr>
```

## Testing with docker and Visual Studio Code

To replicate the environment used in development, you need

1.  [Visual Studio Code](https://code.visualstudio.com/) with the
    [Remote-Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
    extension.
2.  [Docker
    Desktop](https://www.docker.com/products/docker-desktop/alternatives)
    or just docker and docker-compose.

After installing the requisites above,

1.  Clone this repository
    `git clone https://https://github.com/jlaria/glasp.git`
2.  In vscode, go to `File/Open Folder` and open the root folder of this
    project `glasp`.
3.  Then, using the Remote-Containers extension, Reopen in Container (or
    ctrl+p and type `>Remote-Containers:Reopen in Container`)
4.  Wait for everything to install and then go to
    `http://localhost:8787` username:rstudio, password:rstudio.

That’s it! You can use an isolated rstudio-server for development.

See `inst/rstudio` for more details about the environment.
