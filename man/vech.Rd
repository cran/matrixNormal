% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vech.R
\name{vech}
\alias{vech}
\title{Half-Vectorization of a matrix}
\usage{
vech(A, use.Names = TRUE, tol = .Machine$double.eps^0.5)
}
\arguments{
\item{A}{A matrix with m rows and n columns.}

\item{use.Names}{Logical. If TRUE, the names of A are taken to be names of the stacked matrix. Default: TRUE.}

\item{tol}{A numeric tolerance level used to check if a matrix is symmetric. That is, a matrix is symmetric if the difference between the matrix and its transpose is between -\code{tol} and \code{tol}.}
}
\value{
A vector with n(n+1)/2 elements.
}
\description{
Stacks elements of the lower triangle of a numeric symmetric matrix \emph{A}.
}
\details{
For a symmetric matrix \emph{A}, the vectorization of \emph{A} contains more information than necessary. The half-vectorization, denoted \code{vech()}, of a symmetric square n by n matrix \emph{A} is the vectorization of the lower triangular portion.
}
\note{
Unlike other \code{vech()} functions available on CRAN, matrixNormal version may inherit names from matrices to their vectorized forms.
}
\examples{
x <- matrix(c(1, 2, 2, 4),
  nrow = 2, byrow = TRUE,
  dimnames = list(1:2, c("Sex", "Smoker"))
)
print(x)

# Example 1
vech(x)
# If you just want the vectorized form
vech(x, use.Names = FALSE)

# Example 2: If one has NA's
x[1, 2] <- x[2, 1] <- NA
vech(x)
}
