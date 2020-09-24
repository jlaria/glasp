
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Group Linear Algorithm with Sparse Principal decomposition <img src="man/figures/glasp.png" align="right" width="150" />

> A variable selection and clustering method for generalized linear
> models

[![Travis build
status](https://travis-ci.org/jlaria/glasp.svg?branch=master)](https://travis-ci.org/jlaria/glasp)
[![](https://img.shields.io/github/languages/code-size/jlaria/glasp.svg)](https://github.com/jlaria/glasp)
[![](https://img.shields.io/github/last-commit/jlaria/glasp.svg)](https://github.com/jlaria/glasp/commits/master)
[![](https://img.shields.io/badge/devel%20version-0.0.1.9000-blue.svg)](https://github.com/jlaria/glasp)

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
#> Set the global option `yardstick.event_first` to `FALSE` to change this.
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
model <- glasp_model(l1=0.01, l2=0.001, frob=0.0003) %>%
    set_mode("regression") %>%
    set_engine("glasp") %>%
    fit(y~., data)
```

The coefficients can be accessed through the `parsnip` model object
`model$fit$beta`.

``` r
print(model)
#> parsnip model object
#> 
#> Fit time:  57ms 
#> <linear_regression> 
#> $beta
#>           X1           X2           X3           X4           X5           X6 
#>  0.004520053  0.402060668 -0.032605834 -0.369015410 -0.689534760  0.000000000 
#>           X7           X8           X9          X10 
#>  0.000000000  0.023018621  0.000000000  0.000000000 
#> 
#> $intercept
#> [1] -0.1668445
#> 
#> $clusters
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   1   1   1   1   0   1   1   1   1   1 
#> 
#> $info
#> $info$l1
#> [1] 0.01
#> 
#> $info$l2
#> [1] 0.001
#> 
#> $info$frob
#> [1] 3e-04
#> 
#> $info$num_comp
#> [1] 1
```

Now we generate some out-of-sample data to check how the model predicts.

``` r
set.seed(1)
new_data <- simulate_dummy_linear_data()

pred <- predict(model, new_data)
rmse <- rmse_vec(new_data$y, pred$.pred)
```

We obtain a 1.518016 root mean square error.

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

model <- glasp_model(l1 = tune(), 
                     l2 = tune(),
                     frob = tune(), 
                     num_comp = tune()) %>%
         set_mode("classification") %>% 
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
#> ── Iteration 1 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6096 (@iter 0)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=1.11e-06, l2=8.37, frob=3.84e-06, num_comp=8
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5976 (+/-0.0404)
#> 
#> ── Iteration 2 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6096 (@iter 0)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=9.75, l2=6.7, frob=3.72e-06, num_comp=11
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5976 (+/-0.0404)
#> 
#> ── Iteration 3 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6096 (@iter 0)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=0.000441, l2=1.08e-05, frob=0.000219, num_comp=11
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.6112 (+/-0.0586)
#> 
#> ── Iteration 4 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6112 (@iter 3)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=0.0474, l2=1.03e-06, frob=0.00435, num_comp=10
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.6176 (+/-0.0505)
#> 
#> ── Iteration 5 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6176 (@iter 4)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=8.36, l2=1.46e-06, frob=0.0192, num_comp=3
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5976 (+/-0.0404)
#> 
#> ── Iteration 6 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6176 (@iter 4)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=5.25e-06, l2=3.66e-06, frob=6.4, num_comp=7
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5968 (+/-0.0399)
#> 
#> ── Iteration 7 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6176 (@iter 4)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=0.00236, l2=0.00667, frob=0.000906, num_comp=15
#> i Estimating performance
#> ✓ Estimating performance
#> ♥ Newest results:    roc_auc=0.6233 (+/-0.0578)
#> 
#> ── Iteration 8 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6233 (@iter 7)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=8.3e-05, l2=0.000221, frob=0.142, num_comp=15
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.6056 (+/-0.0447)
#> 
#> ── Iteration 9 ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6233 (@iter 7)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=0.128, l2=0.361, frob=4.36e-06, num_comp=2
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5976 (+/-0.0404)
#> 
#> ── Iteration 10 ──────────────────────────────────────────────────────────────────────────────────────────────────────────
#> 
#> i Current best:      roc_auc=0.6233 (@iter 7)
#> i Gaussian process model
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> ✓ Gaussian process model
#> i Generating 5000 candidates
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> * NA -> ...3
#> * NA -> ...4
#> i Predicted candidates
#> i l1=9.48, l2=3.97e-06, frob=9.52, num_comp=15
#> i Estimating performance
#> ✓ Estimating performance
#> ⓧ Newest results:    roc_auc=0.5976 (+/-0.0404)

show_best(hist, metric = "roc_auc")
#> # A tibble: 5 x 10
#>         l1      l2    frob num_comp .iter .metric .estimator  mean     n std_err
#>      <dbl>   <dbl>   <dbl>    <int> <dbl> <chr>   <chr>      <dbl> <int>   <dbl>
#> 1  2.36e-3 6.67e-3 9.06e-4       15     7 roc_auc binary     0.623     4  0.0578
#> 2  4.74e-2 1.03e-6 4.35e-3       10     4 roc_auc binary     0.618     4  0.0505
#> 3  4.41e-4 1.08e-5 2.19e-4       11     3 roc_auc binary     0.611     4  0.0586
#> 4  3.49e-4 5.49e-4 2.68e-5        7     0 roc_auc binary     0.610     4  0.0573
#> 5  8.30e-5 2.21e-4 1.42e-1       15     8 roc_auc binary     0.606     4  0.0447
```

## Install in docker

There is a docker image available with R 3.6.3 and all the dependencies
of `glasp`.

    docker pull jlaria/glasp:0.0.1
    docker run -it jlaria/glasp:0.0.1

Optionally, it can be built using the following `Dockerfile`

    FROM r-base:3.6.3
    LABEL maintainer juank.laria@gmail.com
    
    RUN apt-get update && apt-get install r-cran-dplyr r-cran-rcpparmadillo r-cran-mass r-cran-devtools r-cran-devtools r-cran-ggpubr -y
    
    RUN Rscript -e 'devtools::install_github("jlaria/glasp")'
    
    CMD ["/bin/bash"]
