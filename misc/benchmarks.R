library(bootR2)
## simulate data
set.seed(1)
n <- 500  # number of observations (e.g. cases/sites)
p <- 50  # number of explanatory variables
X <- cbind(1, matrix(rnorm(n * p), n, p))  # model matrix
y <- X %*% rnorm(p + 1, sd = 0.1) + rnorm(n)  # response variable

## compute bootstrap sample and give its summary
yBootR2 <- bootR2(X, y, nBoot = 10000)
summary(yBootR2)

bootR2pureR1 <- function(X, y) {
    n <- nrow(X)
    prmt <- sample.int(n, n, replace = TRUE)
    prmv <- sample.int(n, n, replace = TRUE)
    Xt <- X[prmt, ]
    Xv <- X[prmv, ]
    yt <- y[prmt]
    yv <- y[prmv]
    ## betaHat <- bootR2:::betaHat(Xt, yt)
    betaHat <- qr.coef(qr(Xt), yt)
    fitHat <- Xv %*% betaHat
    SSerr <- sum((yv - fitHat)^2)
    SStot <- sum((yv - mean(yv))^2)
    1 - (SSerr/SStot)
}
bootR2pureR <- function(X, y, nBoot = 1) replicate(nBoot, bootR2pureR1(X, y))
yBootR2pureR <- bootR2pureR(X, y, 10000)
hist(yBootR2pureR, 50, main = "", xlab = expression(paste("Predictive ", R^2)))
hist(yBootR2, 50, main = "", xlab = expression(paste("Predictive ", R^2)))
library(rbenchmark)

benchmark(bootR2(X, y, nBoot = 10000), bootR2pureR(X, y, nBoot = 10000), order = "relative", 
          replications = 10)
