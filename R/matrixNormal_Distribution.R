#'The Matrix Normal Distribution

#'@family distribution
#'@keywords distribution

#'@name matrixNormal_Distribution
NULL
#'@description The density (dmatnorm), cumulative distribution function (CDF, pmatnorm), and generation of 1 random number from the matrix normal (rmatnorm) is produced:
#' A ~ MatNorm_n,p(M, U, V)
#'
#'@details
#'Ideally, both scale matrices are positive-definite. However, they may not appear to be symmetric; you may want to increase the tolerance.
#'
#' These functions rely heavily on this property of matrix normal distribution. Let function koch() refer to the Kronecker product of a matrix. For a n x p matrix A, if \eqn{A ~ MatNorm(M, U, V)}, then \deqn{ vec(A) ~ MVN_{np} (M, Sigma = koch(U,V) ) .} Thus, we can find the probability that Lower < A < Upper by finding the CDF of vec(A), which is given in \code{\link[mvtnorm]{pmvnorm}} function in \pkg{mvtnorm}. Also, we can simulate 1 random matrix A from a matrix normal by sampling vec(A) from \code{\link[mvtnorm]{rmvnorm}} form\pkg{mvtnorm}.

#'@references
#' Iranmanesh, Anis, M. Arashi, and S. M. M. Tabatabaey  On Conditional Applications of Matrix Variate Normal Distribution. \emph{Iranian Journal of Mathematical Sciences and Informatics} 5, no. 2. (November 1, 2010): 33-43. < https://doi.org/10.7508/ijmsi.2010.02.004 >

#'@param A  The numeric n x p matrix that follows the matrix-normal.
#'@param M  The mean n x p matrix that is numeric and real.
#'@param U  The individual scale n x n real positive-definite matrix
#'@param V  The parameter scale p x p  real positive-semidefinite matrix
#'@inheritParams is.symmetric.matrix
            #tol
#'@param use.log Logical; if TRUE, densities d are given as log(d).

#'@examples
#' #Data Used
#' A <- CO2[1:10, 4:5]
#' M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11) )
#' V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
#' V  #Right covariance matrix (2 x 2), say the covariance between parameters.
#' U <- I(10) #Block of left-covariance matrix ( 84 x 84), say the covariance between subjects.
#' #PDF
#' dmatnorm(A, M, U, V )
#' dmatnorm(A, M, U, V, use.log = FALSE)
#'
#' #Generating Probability Lower and Upper Bounds (They're matrices )
#' Lower <- matrix( rep(-1, 20), ncol = 2)
#' Upper <- matrix( rep(3, 20), ncol = 2)
#' Lower; Upper
#' #The probablity that a randomly chosen matrix A is between Lower and Upper
#' pmatnorm( Lower, Upper, M, U, V)
#' #CDF
#' pmatnorm( Lower = -Inf, Upper, M, U, V)
#' #entire domain = 1
#' pmatnorm( Lower = -Inf, Upper = Inf, M, U, V)
#'
#' #Random generation
#' M <- cbind(rnorm(3, 435, 296), rnorm(3, 27, 11) )
#' U <- diag(1, 3)
#' V <- matrix(c(10, 5 ,5, 3), nrow = 2)
#' set.seed(123)
#' rmatnorm(M, U, V)

#' \dontrun{  #M has a different sample size than U; will return an error.
#' M <- cbind(rnorm(4, 435, 296), rnorm(4, 27, 11) )
#' rmatnorm(M, U, V)
#' }

#'@rdname matrixNormal_Distribution
#'@importFrom mvtnorm pmvnorm
#'@export dmatnorm
dmatnorm <- function(A, M, U, V, tol = .Machine$double.eps^0.5, use.log = TRUE){
  n <- nrow(A)
  p <- ncol(A)

  #Checks
  if(is.data.frame(A) ) { A <- as.matrix(A) }
  if( sum(dim(A) == dim(M)) != 2 ) {stop("M must have same dimensions as A.")}
  check_matnorm(M, U, V, tol)

  #The Log Density
  log.dens <- (-n*p/2)*log(2*pi) - p/2*log( det(U)) - n/2*log( det(V) ) +
    -1/2 * tr( solve(U) %*% (A - M) %*% solve(V) %*% t(A - M) )

  #Return
  if( use.log) { return (log.dens) } else { return( exp(log.dens) ) }
}

#'@rdname matrixNormal_Distribution
#'@param Lower	 The n x p matrix of lower limits for CDF
#'@param Upper	 The n x p matrix of upper limits for CDF
#'@inheritParams mvtnorm::pmvnorm
#'@export pmatnorm
pmatnorm <- function( Lower = -Inf, Upper = Inf, M, U, V, tol = .Machine$double.eps^0.5, algorithm  = mvtnorm::GenzBretz(), ...){
  n <- nrow(M)
  p <- ncol(M)

  #Checks
  check_matnorm(M, U, V, tol)

  #Convert the matrices to lower
  if( is.matrix(Lower) ) {
    lower <- vec(Lower)
  } else {
    if( is.vector(Lower) & Lower == -Inf  ) {
      lower <- -Inf
    } else {
      stop( "The lower limit must be a numeric matrix or -Inf.")
    }
  }

  if( is.matrix(Upper) ) {
    upper <- vec(Upper)
  } else {
    if( is.vector(Upper) & Upper == Inf ) {
      upper <- Inf
    } else {
      stop( "The upper limit must be a numeric matrix or Inf.")
    }
  }
  #Calculating the probablity
  prob <- mvtnorm::pmvnorm(lower, upper, mean = vec(M), corr = NULL, sigma = kronecker(U,V),
                           algorithm = algorithm, ...)
  warning("The covariance matrix is standardized. ")

  return(prob)
}

#'@rdname matrixNormal_Distribution
#'@inheritParams mvtnorm::rmvnorm
#'@import mvtnorm
#'@export rmatnorm
rmatnorm <- function(M, U, V, tol = .Machine$double.eps^0.5, method = "chol", pre0.9_9994 = FALSE){
  n <- nrow(M)
  p <- ncol(M)

  #Checks
  check_matnorm( M, U, V, tol)

  #Vectorizing and sampling from rmvnorm
  Sigma <- kronecker(U,V)
  vec.X <- mvtnorm::rmvnorm(1, vec(M), Sigma, method = method, pre0.9_9994)
  X <- matrix(vec.X, nrow = n, ncol = p)
  return(X)
}

# cov(vec(A))  #should be 1


#Check to make sure the parameters in MatrixNormal match.
check_matnorm <- function( M, U, V, tol){
  if(nrow(M) != nrow(U)){stop( "The mean matrix M has different sample size than scale sample size matrix U.")}
  if(ncol(M) != nrow(V)){ stop( "The mean matrix M has different parameters than scale parameter matrix V.")}
  if(!is.positive.definite(U, tol) ){
    stop( "U is not positive definite. Calculation may not be accurate. Possibly raise tolerance"   )
  }
  if(!is.positive.definite(V, tol) ){
   stop( "V is not positive definite. Calculation may not be accurate. Possibly raise tolerance.")
  }
}
