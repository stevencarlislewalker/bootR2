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

```
library(bootR2)
>                                         ## simulate data
> set.seed(1)
> n <- 100 # number of observations (e.g. cases/sites)
> p <- 10  # number of explanatory variables
> X <- cbind(1, matrix(rnorm(n*p), n, p))  # model matrix
> y <- X%*%rnorm(p+1, sd = 0.1) + rnorm(n) # response variable
> 
>                                         ## compute bootstrap sample
>                                         ## and give its summary
> yBootR2 <- bootR2(X, y, nBoot = 10000)
> summary(yBootR2)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.54930  0.02755  0.11130  0.09918  0.18770  0.45510 
```
Note the bootstrap samples with negative R-squared statistics.  This
makes sense for out-of-sample predictions, and indicates that the
model predictions are more variable than the validation data
themselves.


Limitations
-----------

1. currently only least-squares fits are permitted 
2. currently only sum-of-squares-based R-squared statistics are computed
