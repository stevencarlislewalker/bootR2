library(bootR2)
                                        ## simulate data
set.seed(1)
n <- 100 # number of observations (e.g. cases/sites)
p <- 10  # number of explanatory variables
X <- cbind(1, matrix(rnorm(n*p), n, p))  # model matrix
y <- X%*%rnorm(p+1, sd = 0.1) + rnorm(n) # response variable

                                        ## compute bootstrap sample
                                        ## and give its summary
yBootR2 <- bootR2(X, y, nBoot = 10000)
summary(yBootR2)


bR2 <- bootR2samp(cbind(1, X[,4,drop=FALSE]), y, 10000)
bprR2 <- bootR2pred(cbind(1, X[,4,drop=FALSE]), y, 10000)
par(mfrow = c(3, 1))
plot(X[,4], y)
hist(bR2, 50)
hist(bprR2, 50)

m <- lm(y ~ 0 + ., as.data.frame(X))
m0 <- lm(y ~ 0 + X[,4])
(sm <- summary(m0))
R2(X, y)
sm$adj
mean(bprR2)
mean(bprR2 < 0)
mean(bprR2 < sm$adj)



Xv <- cbind(1, rnorm(n))
yv <- Xv%*%c(1, -0.1) + rnorm(n)
betaHat(X, y)
rUnif(n)
bootPerm(n)

order_(y)
order(y)

bootPerm(n)
shuffleMatrix(X, sample(n)-1)
shuffleVector(y, sample(n)-1)
BB <- bootCoef(X, y, 10000)
plot(BB)


hist(bootR2(X, y, 10000))
abline(v = R2(X, y), col = 'red')
abline(v = R2pred(X, y, Xv, yv), col = 'red')


hist(bootR2pred(X, y, 10000))
abline(v = R2(X, y), col = 'red')
abline(v = R2pred(X, y, Xv, yv), col = 'red')



R2(X, y)
R2pred(X, y, Xv, yv)

hist(bootR2(X, y, 10000))
hist(bootR2pred(X, y, 10000))


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
replicate(10000, bootR2pureR1(X, y))

library(vegan)
data(dune)
data(dune.env)
dune.Manure <- rda(dune ~ Manure, dune.env)
dune.mlm <- as.mlm(dune.Manure)
model.response(model.frame(dune.mlm))
model.matrix(dune.mlm)

model.response(model.frame(dune.Manure))

model.frame(formula(dune.Manure$call), )
eval(dune.Manure$call$data)
