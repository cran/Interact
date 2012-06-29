\name{Interact-package}
\alias{Interact-package}
\alias{Interact}
\docType{package}
\title{
Test marginal interactions for a model with binary response
}
\description{
This package tests all potential marginal interactions in a binary response model, and estimates false discovery rates at each potential cutoff.
}
\details{
\tabular{ll}{
Package: \tab Interact\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-29\cr
License: \tab GPL-2\cr
}
Very straightforward to use. Accepts data in the form \code{X,y} and tests all marginal interactions. Only 2 functions:
\code{interact}
\code{plot.interact}
}
\author{Noah Simon, and Robert Tibshirani|

Maintainer: Noah Simon <nsimon@stanford.edu>

}
\references{Simon, N. and Tibshirani, R. (2012) \emph{A Permutation Approach to Testing Marginal Interactions in Many Dimensions}, \url{http://www-stat.stanford.edu/~nsimon/TMIcor.pdf}
}

\keyword{interactions}
\keyword{FDR}


