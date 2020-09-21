SSL Validation
================

## Overview

Phenotyping algorithms based on Electronic Health Records (EHR) aim to
identify a patient’s disease status using the information in the health
record.

Algorithm evaluation in EHR data is often based on a small number of
gold-standard labeled data. Two main issues are the difficulty in
selecting the appropriate cut and the high variance in the accuracy
parameter estimates.

This package provides a semi-supervised learning approach that
incorporates unlabeled data into the estimation of the receiver
operating characteristic (ROC) curve parameters to address these issues.

The data consists of :

  - a small validation set with algorithm score S and label Y, which
    takes values 0 or 1

  - a large unlabeled set containing only the algorithm score S.

The main function `roc.semi.superv` takes as arguments S and Y, where Y
contains a large amount of missing values, and provides a
semi-supervised estimation of the ROC parameters.

The function `roc.superv` takes as arguments S and Y from the small
validation set only and provides a classic supervised method to estimate
the ROC parameters.

## Installation

Install development version from GitHub.

``` r
# install.packages("remotes")
remotes::install_github("celehs/SSL.validation")
```

Load the package into R.

``` r
library(SSL.validation)
```

## Simulated Example

``` r
set.seed(1234)
dat <- read.csv("https://raw.githubusercontent.com/celehs/SSL.validation/master/data-raw/data.csv")
```

``` r
p.0 <- mean(dat$Y , na.rm = TRUE)
id.v <- which(is.na(dat$Y) != 1)
dat.v <- dat[id.v, ] # Labeled Data 
```

### Semi-Supervised Learning (SSL)

``` r
system.time(res.ssl <- roc.semi.superv(dat$S,dat$Y))
```

    ## [1] 263

    ##    user  system elapsed 
    ##  24.285   0.856  29.674

``` r
auc.ssl <- res.ssl$auc
roc.ssl <- res.ssl$roc
auc.ssl
```

    ## [1] 0.7467197

``` r
tail(roc.ssl)
```

    ##          cut  p.pos  fpr       tpr       ppv       npv
    ## [94,] 0.0471 0.9614 0.94 0.9909063 0.2591625 0.9518820
    ## [95,] 0.0386 0.9688 0.95 0.9930972 0.2574391 0.9554319
    ## [96,] 0.0312 0.9785 0.96 0.9949573 0.2559512 0.9597191
    ## [97,] 0.0215 0.9885 0.97 0.9973425 0.2540214 0.9691953
    ## [98,] 0.0115 1.0000 0.98 0.9997318 0.2520541 0.9941873
    ## [99,] 0.0115 1.0000 0.98 0.9997318 0.2520541 0.9941873

### Supervised Learning (SL)

``` r
system.time(res.sl <- roc.superv(dat.v$S,dat.v$Y))
```

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ##    user  system elapsed 
    ##   0.006   0.000   0.006

``` r
auc.sl <- res.sl$auc
roc.sl <- res.sl$roc
colnames(roc.sl) <- colnames(roc.ssl)
auc.sl
```

    ## [1] 0.7016632

``` r
tail(roc.sl)
```

    ##               cut p.pos       fpr tpr       ppv npv
    ## [94,] 0.015460938  0.96 0.9459459   1 0.2708333   1
    ## [95,] 0.014389420  0.96 0.9459459   1 0.2708333   1
    ## [96,] 0.013239677  0.98 0.9729730   1 0.2653061   1
    ## [97,] 0.010720983  0.98 0.9729730   1 0.2653061   1
    ## [98,] 0.008202289  0.98 0.9729730   1 0.2653061   1
    ## [99,] 0.006568542  0.99 0.9864865   1 0.2626263   1
