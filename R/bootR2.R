##' Bootstrap R-squared
##'
##' @param object Object for calculating R-squared
##' @return Vector of R-squared bootstrap samples
##' @export
bootR2 <- function(object, ...) UseMethod("bootR2")

bootR2.default <- function(object, ...) {
    msg <- paste("No bootR2 method for objects of class",
                 class(object)[1])
    stop(msg)
}

bootR2.lm <- function(object, nBoot = 1, predictive = TRUE, ...) {
    X <- model.matrix(object)
    y <- model.response(model.frame(object))
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootFunc(X, y, nBoot)
}

bootR2.matrix <- function(object, y, nBoot = 1, predictive = TRUE, ...) {
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootR2pred(object, y, nBoot)
}

bootR2.numeric <- function(object, X, nBoot = 1, predictive = TRUE, ...) {
    bootFunc <- ifelse(predictive, bootR2pred, bootR2samp)
    bootR2pred(X, object, nBoot)
}
