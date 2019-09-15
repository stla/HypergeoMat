#' Incomplete Beta function of a matrix argument
#'
#' @description Evaluates the incomplete Beta function of a matrix argument.
#'
#' @param m truncation weight, a positive integer
#' @param a,b real or complex parameters with \code{Re(a)>(p-1)/2},
#' \code{Re(b)>(p-1)/2}, where \code{p} is the dimension (the order of the matrix)
#' @param x either a real positive symmetric matrix or a complex positive
#' Hermitian matrix "smaller" than the identity matrix (i.e. \code{I-x} is positive),
#' or a numeric or complex vector, the eigenvalues of the matrix
#'
#' @return A real or a complex number.
#' @export
#'
#' @note The eigenvalues of a real symmetric matrix or a complex Hermitian
#' matrix are always real numbers, and moreover the are positive under the
#' constraints on \code{x}.
#' However we allow to input a numeric or complex vector \code{x}
#' because the definition of the function makes sense for such a \code{x}.
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a scalar x, this is the incomplete Beta function:
#' a <- 2; b <- 3
#' x <- 0.75
#' IncBeta(m = 15, a, b, x)
#' gsl::beta_inc(a, b, x)
IncBeta <- function(m, a, b, x){
  if(is.matrix(x)){
    stopifnot(isSymmetricPositive(x))
    p <- nrow(x)
    stopifnot(isSymmetricPositive(diag(p)-x))
  }else if(is.atomic(x)){
    p <- length(x)
  }else{
    stop("Invalid `x` argument")
  }
  stopifnot(
    is.numeric(a) || is.complex(a),
    Re(a) > (p-1)/2,
    is.numeric(b) || is.complex(b),
    Re(b) > (p-1)/2
  )
  if(is.matrix(x)){
    DET <- det(x)
  }else{
    DET <- prod(x)
  }
  DET^a / mvbeta(a, b, p) / a *
    hypergeomPFQ(m, c(a, -b+(p+1)/2), a+(p+1)/2, x)
}
