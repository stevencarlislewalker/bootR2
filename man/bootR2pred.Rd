\name{bootR2pred}
\alias{bootR2pred}
\title{The main functions for the package}
\usage{
bootR2pred(X, y, nBoot)
}
\arguments{
  \item{X}{model matrix}

  \item{y}{model response}

  \item{nBoot}{number of bootstrap samples}
}
\value{
vector of \code{nBoot} R-squared values
}
\description{
\code{bootR2pred} gives out-of-sample bootstrap samples of
R-squared, whereas \code{bootR2samp} gives within-sample
R-squared bootstrap samples
}

