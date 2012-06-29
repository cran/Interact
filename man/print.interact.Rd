\name{print.interact}
\alias{print.interact}

\title{Prints Basic output for interaction object
}
\description{
This prints the function call and first 10 significant interactions for an object output by the function \code{interact}
}
\usage{
\method{print}{interact}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A variable of class \code{interact} --- the output of function \code{interact}
}
  \item{\dots}{Other arguments which the user would like to pass to the general \code{print} function
  
}
}

\references{Simon, N. and Tibshirani, R. (2012) \emph{A Permutation Approach to Testing Marginal Interactions in Many Dimensions}, \url{http://www-stat.stanford.edu/~nsimon/TMIcor.pdf}
}

\author{Noah Simon, and Robert Tibshirani}

\seealso{ \code{interact}, \code{plot}
}


\keyword{interactions}
\keyword{FDR}
