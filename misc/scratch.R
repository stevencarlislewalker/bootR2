library(bootR2)
n <- 100
p <- 10
X <- matrix(rnorm(n*p), n, p)
y <- X%*%rnorm(p, sd = 0.15) + rnorm(n)
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

