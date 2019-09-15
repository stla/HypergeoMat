#' Type one Bessel function of Herz
#'
#' @description Evaluates the type one Bessel function of Herz.
#'
#' @param m truncation weight, a positive integer
#' @param nu the order parameter, real or complex number with \code{Re(nu)>-1}
#' @param x either a real symmetric matrix, a Hermitian complex matrix,
#' or a numeric or complex vector, the eigenvalues of the matrix
#'
#' @return A real or complex number.
#' @export
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a scalar x, the relation with the Bessel J-function:
#' t <- 2
#' nu <- 3
#' besselJ(t, nu)
#' Bessel(m=15, t^2/4, nu) * (t/2)^nu
#' # it also holds for a complex variable:
#' t <- 1 + 2i
#' Bessel::BesselJ(t, nu)
#' Bessel(m=15, t^2/4, nu) * (t/2)^nu
Bessel <- function(m, x, nu){
  stopifnot(is.numeric(nu) || is.complex(nu), Re(nu) > -1)
  if(is.matrix(x)){
    p <- nrow(x)
  }else{
    p <- length(x)
  }
  hypergeomPFQ(m, NULL, nu+(p+1)/2, -x) / mvgamma(nu+(p+1)/2, p)
}
