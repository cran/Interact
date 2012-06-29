\name{plot.interact}
\alias{plot.interact}

\title{Plots FDR estimates for interact object
}
\description{
This plots an estimated FDR curve for an object output by the function \code{interact}
}
\usage{
\method{plot}{interact}(x, numInteractions = nrow(x$interaction.ordered), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A variable of class \code{interact} --- the output of function \code{interact}
}
  \item{numInteractions}{The number of most significant interactions one would like plotted
}
  \item{\dots}{Other arguments which the user would like to pass to the general \code{plot} function
  
}
}

\references{Simon, N. and Tibshirani, R. (2012) \emph{A Permutation Approach to Testing Marginal Interactions in Many Dimensions}, \url{http://www-stat.stanford.edu/~nsimon/TMIcor.pdf}
}

\author{Noah Simon, and Robert Tibshirani}

\seealso{ \code{interact}, \code{print}
}


\keyword{interactions}
\keyword{FDR}
