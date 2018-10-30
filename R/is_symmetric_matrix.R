#' Is a matrix symmetric or positive-definite?
#'
#' @name is.symmetric.matrix
#' @family statistics
#' @keywords matrix
#'
#' @description Determine if a matrix is square, symmetric, positive-definite, or positive-semi-definite.
#' @details
#' A tolerance is added to indicate if a matrix is approximately symmetric. If the matrix is not symmetric, a message as well as the top of the matrix is printed.
#' \itemize{
#' \item \code{is.symmetric.matrix} returns TRUE if A is a numeric symmetric square matrix and FALSE otherwise. A matrix is symmetric if the difference between A and its transpose is less than \emph{tol}.
#' \item \code{is.positive.semi.definite} returns TRUE if a square symmetric real matrix A is positive semi-definite. A matrix is positive semi-definite if its smallest eigenvalue is greater than or equal to zero.
#' \item \code{is.positive.definite} returns TRUE if a square symmetric real matrix A is positive-definite.A matrix is positive-definite if its smallest eigenvalue is greater than zero. }

#' @note
#' Functions adapted from Frederick Novomestky's \pkg{matrixcalc} package in order to implement \code{rmatnorm} function.  I changed argument x to A to reflect usual matrix notation. For \code{is.symmetric}, I added a tolerance so that A is symmetric even provided small differences between A and its transpose. Useful for rmatnorm function, which was used repeatedly to generate matrixNormal random variates in a Markov chain. For \code{is.positive.semi.definite} and \code{is.positive.definite},  I also saved time by avoiding a for loop and instead calculating the minimum of eigenvalues.

#' @param A Numeric matrix with no missing values.
#' @param tol A numeric tolerance level used to check if a matrix is symmetric; that is if the difference between it and its transpose is between -tol and tol.

# May want to produce a warning instead of stopping.

#' @examples
#' ## Example 0: Not square matrix
#' B <- matrix( c( 1, 2, 3, 4, 5, 6 ), nrow=2, byrow=TRUE )
#' B
#' is.square.matrix(B)
#'
#' ## Example 1: Not a matrix. should get an error.
#' \dontrun{ df <- as.data.frame( matrix( c( 1, 2, 3, 4, 5, 6 ), nrow=2, byrow=TRUE ) )
#' df
#' is.square.matrix(df)
#' }
#'
#' ## Example 2: Not Symmetric & Compare against matrixcalc
#' if( !requireNamespace( "matrixcalc", quietly = TRUE)) { install.packages("matrixcalc") }
#'   F <- matrix( c( 1, 2, 3, 4 ), nrow=2, byrow=TRUE ); F
#' is.square.matrix(F)
#' is.symmetric.matrix( F )   #should be FALSE
#' matrixcalc::is.symmetric.matrix(F)
#' #Another Symmetric Test found in base. Because of this, is.symmetric() may not be needed
#' isSymmetric.matrix(F)
#'
#' ## Example 3: Symmetric but negative-definite. same test of functions
#' ##' eigenvalues are  3 -1
#' G <- matrix( c( 1, 2, 2, 1 ), nrow=2, byrow=TRUE ); G
#' is.symmetric.matrix( G )
#' matrixcalc::is.symmetric.matrix(G)
#' isSymmetric.matrix(G)
#' is.positive.definite(G) #FALSE
#' is.positive.semi.definite(G) #FALSE
#'
#' ## Example 4: positive definite matrix
#' #' eigenvalues are 3.4142136 2.0000000 0.585786
#' Q <- matrix( c( 2, -1, 0, -1, 2, -1, 0, -1, 2 ), nrow=3, byrow=TRUE )
#' is.symmetric.matrix(Q)
#' is.positive.definite( Q )
#'
#' ## Example 5: identity matrix is always positive definite
#' I <- diag( 1, 3 )
#' is.square.matrix(I) #TRUE
#' is.symmetric.matrix(I) #TRUE
#' is.positive.definite( I ) #TRUE

#' @importFrom utils head

#' @rdname is.symmetric.matrix
#' @export is.square.matrix
#Is matrix square? (A must be a matrix. ).
#Adapted from matrixcalc::is.square.matrix. Argument name changed to be consistent.
is.square.matrix <- function (A){
  if (!is.numeric(A))   stop("A is not a numeric matrix")
  if (!is.matrix(A)) stop( "A is not a matrix")
  is.square <- nrow(A) == ncol(A)
  return( is.square )
}

#' @rdname is.symmetric.matrix
#' @export is.symmetric.matrix
# Is matrix is symmetric? A must be a numeric square matrix with no missing values.
# ??? Prefer the error to call name of matrix.
is.symmetric.matrix <- function (A, tol = .Machine$double.eps^0.5){
  #Checks: A must be a numeric square matrix with no missing values.
  if(anyNA(A)){  return(NA)}
  stopifnot (is.square.matrix(A) ) #stop( "A is not a square matrix")

  #Edited from matrixcalc (PH)
  #Is A and t(A) equal within a tolerance?
  total.abs <-  sum ( abs(A - t(A) ) )
   if(  total.abs   < tol ) {    # to avoid being exactly equal
    okay <- TRUE
    #print("sample Variance Matrix is symmetric")
  } else {
    print("A is not symmetric. Top of the matrix: ")
    print(utils::head(A))
    okay <- FALSE
  }
  # cat("sum( abs(A - t(A)) : ", total.abs, "\n")
  # cat("Total Absolute Difference between Matrix & It's Transpose: ", total.abs, "\n")
  return(okay)
}

#Find eigenvalues of a symmetric matrix A with no missing values.
#Paul Hargarten added this function.
find.eval <- function(A, tol = .Machine$double.eps^0.5){
  #Check if A is symmetric.
  is.symm <- is.symmetric.matrix(A, tol)
  if(is.symm){
    #If A is symmetric, find eigenvalues.
    eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    #Adjust small eigenvalues to be 0  #(Edited from for loop)
    eigenvalues <- ifelse( abs(eigenvalues) < tol, 0, eigenvalues)
  } else{
    #PAUL HARAGARTEN CHANGED HERE TO A WARNING.  ???????
    stop( "No eigenvalues can be imputed; A is not symmetric.")
  }
  return(eigenvalues)
}

#' @rdname is.symmetric.matrix
#' @export
is.positive.semi.definite <- function (A, tol = .Machine$double.eps^0.5){
  #Positive semi-definite matrix have non-negative eigenvalues.
  eigenvalues <- find.eval(A, tol)
  pos.semi <- if ( min(eigenvalues) >= 0 ) {  TRUE } else { FALSE }
  return(pos.semi)
}

#' @rdname is.symmetric.matrix
#' @export
is.positive.definite <- function (A, tol = .Machine$double.eps^0.5){
  #Positive definite matrix have positive e-values.
  eigenvalues <- find.eval(A, tol)
  pos <- if ( min(eigenvalues) > 0 ) {  TRUE } else { FALSE }
  return(pos)
}


#  #Checks
#USE THIS INSTEAD ????
#Edited from tmvtnorm so that the symmetric test includes eigenvalues instead of the deteriminant.
checkSymmetricPositiveDefinite <- function (x, name = "sigma"){
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  if (any(diag(x) <= 0)) {
    stop(sprintf("%s all diagonal elements must be positive",
                 name))
  }
  min.eval <- min(eigen(x)$values)    #edited from tmvnorm function
  if (det(x) <= 0 & min.eval(x) <0) {
    stop(sprintf("%s must be positive definite", name))
  }
}
