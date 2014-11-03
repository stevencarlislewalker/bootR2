library(bootR2)

N <- 100
n <- 20
m <- 10
p <- 10

Y <- matrix(rpois(N*m, 1), N, m)
X <- matrix(rnorm(N*p), N, p)
iterSimExperiment(X, Y, 5, 5, 10, N)

prm <- bootPerm(n, N)
shuffleMatrix(X, prm)

bootR2(X, Y)


R2pred(, X,
       hellinger(shuffleMatrix(as.double(Y), prm)),
       hellinger(as.double(Y)))

shuffleMatrix(X, 0:10)



simExperiment(shuffleMatrix(X, prm), X,
              hellinger(as.double(shuffleMatrix(as.double(Y), prm))),
              hellinger(as.double(Y)), 2, 2)




bootR2(X, hellinger(Y), nBoot = 10)

storage.mode(X) <- "double"
sqrt(sweep(X, 1, rowSums(X), "/"))
decostand(X, method = "hellinger")
hellinger(X)

