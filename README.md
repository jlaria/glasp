
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
#> Fit time:  33ms 
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
#> i l1=5.9e-06, l2=1.84e-06, frob=9.45, num_comp=7
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.5587 (+/-0.0599)
#> 
#> ── Iteration 2 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5587 (@iter 1)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=7.96e-05, l2=1.07e-06, frob=0.00269, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.5618 (+/-0.0784)
#> 
#> ── Iteration 3 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.0089, l2=1e-06, frob=0.000831, num_comp=6
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5598 (+/-0.085)
#> 
#> ── Iteration 4 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=9.95, l2=0.572, frob=1.53, num_comp=10
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 5 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.105, l2=0.511, frob=3.18e-06, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 6 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.000432, l2=5.01, frob=2.29, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 7 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=0.614, l2=1.59, frob=9.66, num_comp=6
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 8 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=9.03, l2=0.0107, frob=2.19e-06, num_comp=8
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 9 ─────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=7.64, l2=2.4e-05, frob=1.67, num_comp=17
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5579 (+/-0.0608)
#> 
#> ── Iteration 10 ────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.5618 (@iter 2)
#> i Gaussian process model
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> i Predicted candidates
#> i l1=4.62e-05, l2=1.21e-05, frob=0.00479, num_comp=4
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5541 (+/-0.073)

show_best(hist, metric = "roc_auc")
#> # A tibble: 5 x 11
#>           l1        l2      frob num_comp .metric .estimator  mean     n std_err
#>        <dbl>     <dbl>     <dbl>    <int> <chr>   <chr>      <dbl> <int>   <dbl>
#> 1 0.0000796    1.07e-6   2.69e-3        4 roc_auc binary     0.562     4  0.0784
#> 2 0.00890      1.00e-6   8.31e-4        6 roc_auc binary     0.560     4  0.0850
#> 3 0.00000590   1.84e-6   9.45e+0        7 roc_auc binary     0.559     4  0.0599
#> 4 0.0000248    6.69e-6   1.22e-1       17 roc_auc binary     0.558     4  0.0683
#> 5 0.376        1.59e-2   8.53e-6       13 roc_auc binary     0.558     4  0.0608
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
#>           X1           X2           X3           X4           X5           X6 
#>  0.321530396  0.414310760 -0.007403158  0.088009349  1.765738321  0.007941341 
#>           X7           X8           X9          X10 
#> -0.057536219  0.107952601  0.191754532  0.000000000 
#> 
#> $model$intercept
#> [1] 0.07803997
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
#> [1] 2.429158e-01 9.376927e-02 2.219503e-02 6.730942e-03 2.312142e-03
#> [6] 8.936364e-04 3.478873e-04 1.338193e-04 6.008179e-05
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
#> # A tibble: 5 x 10
#>           l1        l2      frob num_comp .metric .estimator  mean     n std_err
#>        <dbl>     <dbl>     <dbl>    <int> <chr>   <chr>      <dbl> <int>   <dbl>
#> 1 0.000188   0.0000443   1.09e-6        2 roc_auc binary     0.864     4  0.0185
#> 2 0.00129    0.0000167   7.15e-6        5 roc_auc binary     0.862     4  0.0183
#> 3 0.0000837  0.0341      5.07e-5        7 roc_auc binary     0.843     4  0.0214
#> 4 0.00000367 0.000590    2.86e+0       20 roc_auc binary     0.831     4  0.0388
#> 5 0.0000138  0.00197     3.46e-1       10 roc_auc binary     0.826     4  0.0368
#> # … with 1 more variable: .config <chr>
```

## Install in docker

See `inst/rstudio`
