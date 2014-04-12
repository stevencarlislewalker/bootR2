bootR2
======

`R` package for computing fast `Rcpp`-based predictive R-squares, via the nonparametric bootstrap.  Heavy use is made of `RcppEigen`.

Example
-------

```
> library(bootR2)
> set.seed(1)
> n <- 100
> p <- 10
> X <- cbind(1, matrix(rnorm(n*p), n, p))
> y <- X%*%rnorm(p+1, sd = 0.1) + rnorm(n)
> yBootR2 <- bootR2(X, y, nBoot = 10000)
> summary(yBootR2)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.54930  0.02755  0.11130  0.09918  0.18770  0.45510 
```
