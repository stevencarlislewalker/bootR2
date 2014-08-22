##' Bootstrap R-squared
##'
##' A package for computing fast Rcpp-based predictive R-squares, via
##' the nonparametric bootstrap.  Heavy use is made of RcppEigen.
##' Predictive R-squared statistics use two different samples, one for
##' training a linear model and one for validating it.  The
##' \code{bootR2} package draws two nonparametric bootstrap samples to
##' be used as training and validation data for calculating a
##' predictive R-squared. Such a bootstrap distribution provides a
##' better assessment of the proportion of variation that the fitted
##' model could explain in other samples from the population.  NOTE:
##' the adjusted R-squared does *not* provide this assessment.
##'
##' @param object Object for calculating R-squared
##' @param ... Not currently used
##' @return Vector of R-squared bootstrap samples
##' @rdname bootR2
##' @aliases bootR2-package
##' @export
bootR2 <- function(object, ...) UseMethod("bootR2")

##' @rdname bootR2
##' @export
bootR2.default <- function(object, ...) {
    msg <- paste("No bootR2 method for objects of class",
                 class(object)[1])
    stop(msg)
}

##' @param nBoot Number of bootstrap samples
##' @param predictive Should the predictive (i.e. double) bootstrap be
##' used?  If \code{TRUE}, each of the \code{nBoot} bootstrap samples
##' actually consists of two samples: one for training the linear
##' model and one for validating it.  Therefore, in this case, the
##' bootstrap distribution may be negative when out-of-sample
##' predictions generate too much variance relative to the bias they
##' avoid.  If \code{FALSE}, the same sample is used for both training
##' and validation.
##' @rdname bootR2
##' @export
bootR2.lm <- function(object, nBoot = 1, predictive = TRUE, ...) {
    X <- model.matrix(object)
    y <- model.response(model.frame(object))
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootFunc(X, y, nBoot)
}

##' @rdname bootR2
##' @export
bootR2.matrix <- function(object, y, nBoot = 1, predictive = TRUE, ...) {
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootFunc(object, y, nBoot)
}

##' @rdname bootR2
##' @export
bootR2.numeric <- function(object, X, nBoot = 1, predictive = TRUE, ...) {
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootFunc(X, object, nBoot)
}

##' Adjusted R-square statistics
##'
##' Various adjustments to sample R-square statistics.  Currently only
##' adjusted and generalised cross validation adjustments are
##' available.  TODO:  add final prediction error
##'
##' @param R2 a vector of R-square statistics
##' @param n length of the response
##' @param p number of columns of the response matrix
##' @return an adjusted version of \code{R2}
##' @rdname adjusted
##' @export
adj <- function(R2, n, p) 1 - ((n - 1) / (n - p)) * (1 - R2)

##' @rdname adjusted
##' @export
gcv <- function(R2, n, p) 1 - (n / (n - p)) * (1 - adj(R2, n, p))

