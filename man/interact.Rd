\name{interact}
\alias{interact}
\title{Test marginal interactions for a model with binary response
}
\description{
Tests all potential marginal interactions, and estimates false discovery rates at each potential cutoff
}
\usage{
interact(x, y, z = NULL, numPerm = 100, numFDR = 1000, method = "Pearson", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An n-by-p matrix of covariates - observations in rows, features in columns.
}
  \item{y}{
An n-vector of class labels taking on two values (eg 0,1 or A,B)
}
  \item{z}{
An optional secondary n-by-q matrix of covariates - observations in
rows, features in columns.
}
  \item{numPerm}{
The number of permutations to run
}
  \item{numFDR}{
The number of marginal interactions you would like to estimate FDR for
--- default is 1000 (more interactions can increase runtime).
}
  \item{method}{
A string, either "Pearson" or "Spearman", indicating which type of correlation is to be used.
}
  \item{verbose}{
A boolean flag indicating whether current permutation number should be output.
}
}
\details{
A correlation matrix is constructed for each class (according to
method). The function then apply a fisher transformation to these values and takes
their difference. These values are ordered, and permutations are used to
assess false discovery rate estimates. If no Z matrix is included then
all pairwise correlations are considered for variables in X. If a Z
matrix is included then only correlations between X and Z variables are considered.
}
\value{
\item{interaction.ordered }{A dataframe of the \code{numFDR} most significant marginal interactions (ordered from most significant to least significant). The first two columns indication the interaction and the third column gives an estimated q-value (False Discovery Rate).
}
\item{internals }{Variables used interally for methods relating to interact}
}
\references{Simon, N. and Tibshirani, R. (2012) \emph{A Permutation Approach to Testing Marginal Interactions in Many Dimensions}, \url{http://www-stat.stanford.edu/~nsimon/TMIcor.pdf}
}

\author{Noah Simon, and Robert Tibshirani}

\seealso{ \code{plot}, \code{print}
}
\examples{
set.seed(5)

n <- 100
p <- 10
s <- 5

X1 <- cbind(matrix(rnorm(n*s), ncol = s) + rnorm(n), matrix(rnorm(n*(p-s)), ncol = (p-s)))
X2 <- matrix(rnorm(n * p), ncol = p)

X <- rbind(X1, X2)

colnames(X) <- c("a","b","c","d","e","f","g","h","i","j")
y <- c(rep("y",n),rep("n",n))

fit <- interact(X,y)
print(fit)
plot(fit)

## Bigger Example (restricting the number of top interactions to consider)
## Not run:
\dontrun{
n <- 300
p <- 500
s <- 10

X1 <- cbind(matrix(rnorm(n*s), ncol = s) + rnorm(n), matrix(rnorm(n*(p-s)), ncol = (p-s)))
X2 <- matrix(rnorm(n * p), ncol = p)
X <- rbind(X1, X2)

y <- c(rep("y",n),rep("n",n))

fit <- interact(X,y, numFDR = 50)  
## Restricts the number of most significant interactions to consider to 50
print(fit)
plot(fit)
}

## Example Comparing (simulated) Genes and Enviromental Variables

\dontrun{
n <- 100
p1 <- 100
p2 <- 10


Genes <- matrix(rnorm(n * p1), ncol = p1)
Environment <- matrix(rnorm(n * p2), ncol = p2)

y <- c(rep("y",n/2),rep("n",n/2))

fit <- interact(X = Genes,y, Z = Environment, numFDR = 50)  
## Restricts the number of most significant interactions to consider to 50
print(fit)
plot(fit)
}
## End(**Not run**)


}
\keyword{interactions}
\keyword{FDR}
