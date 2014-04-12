bootR2
======

`R` package for computing fast `Rcpp`-based predictive R-squares, via
the nonparametric bootstrap.  Heavy use is made of `RcppEigen`.

Predictive R-squared statistics use two different samples, one for
training a linear model and one for validating it.  The `bootR2`
package draws two nonparametric bootstrap samples to be used as
training and validation data for calculating a predictive R-squared.
Such a bootstrap distribution provides a better assessment of the
proportion of variation that the fitted model could explain in other
samples from the population.  NOTE: the adjusted R-squared does *not*
provide this assessment.


Why?
----

The technometrics literature has argued for quite a while that the
adjusted R-squared statistic is biased as an estimate of the ability
of a fitted model to explain variation in the population.  However,
constructing an analytical bias correction is: 

1. challenging
2. has resulted in a large number of proposals, none of which dominate the
others in all settings
3. do not readily yield estimates of uncertainty in the point estimates

The nonparametric bootstrap is a good alternative, but obviously
slower.  By using `RcppEigen`, we are able to compute fast approximate
sampling distributions of out-of-sample predictive R-squared
statistics.


Example
-------


```r
library(bootR2)
## simulate data
set.seed(1)
n <- 100  # number of observations (e.g. cases/sites)
p <- 10  # number of explanatory variables
X <- cbind(1, matrix(rnorm(n * p), n, p))  # model matrix
y <- X %*% rnorm(p + 1, sd = 0.1) + rnorm(n)  # response variable

## compute bootstrap sample and give its summary
yBootR2 <- bootR2(X, y, nBoot = 10000)
summary(yBootR2)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -0.5490  0.0276  0.1110  0.0992  0.1880  0.4550
```

```r
hist(yBootR2, 50, main = "", xlab = expression(paste("Predictive ", R^2)))
```

![plot of chunk README_example](figure/README_example.png) 

Note the bootstrap samples with negative R-squared statistics.  This
makes sense for out-of-sample predictions, and indicates that the
model predictions are more variable than the validation data
themselves.


```r
bootR2pureR1 <- function(X, y) {
    n <- nrow(X)
    prmt <- sample.int(n, n, replace = TRUE)
    prmv <- sample.int(n, n, replace = TRUE)
    Xt <- X[prmt, ]
    Xv <- X[prmv, ]
    yt <- y[prmt]
    yv <- y[prmv]
    betaHat <- bootR2:::betaHat(Xt, yt)
    fitHat <- Xv %*% betaHat
    SSerr <- sum((yv - fitHat)^2)
    SStot <- sum((yv - mean(yv))^2)
    1 - (SSerr/SStot)
}
bootR2pureR <- function(X, y, nBoot = 1) replicate(nBoot, bootR2pureR1(X, y))
yBootR2pureR <- bootR2pureR(X, y, 10000)
hist(yBootR2pureR, 50, main = "", xlab = expression(paste("Predictive ", R^2)))
```

![plot of chunk pureR](figure/pureR.png) 



```r
library(rbenchmark)
benchmark(bootR2(X, y, nBoot = 10000), bootR2pureR(X, y, nBoot = 10000), order = "relative", 
    replications = 10)
```

```
##                               test replications elapsed relative user.self
## 1      bootR2(X, y, nBoot = 10000)           10   2.228    1.000     2.174
## 2 bootR2pureR(X, y, nBoot = 10000)           10   7.468    3.352     7.459
##   sys.self user.child sys.child
## 1    0.054          0         0
## 2    0.009          0         0
```


Limitations
-----------

1. currently only least-squares fits are permitted 
2. currently only sum-of-squares-based R-squared statistics are computed
3. currently only nonparametric bootstraps are allowed

Interested in helping remove these limitations?
