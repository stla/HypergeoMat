#' Type one Bessel function of Herz
#'
#' @description Evaluates the type one Bessel of Herz.
#'
#' @param m truncation weight, a positive integer
#' @param gamma the gamma parameter, real or complex number with \code{Re(gamma)>-1}
#' @param x either a real symmetric matrix, a Hermitian complex matrix,
#' or a numeric or complex vector, the eigen values of the matrix
#'
#' @return A real or complex number.
#' @export
#'
#' @references A. K. Gupta and D. K. Nagar.
#' \emph{Matrix variate distributions}. Chapman and Hall, 1999.
#'
#' @examples # for a single x, the relation with the Bessel J-function:
#' t <- 2
#' gamma <- 3
#' besselJ(t, gamma)
#' Bessel(m=15, t^2/4, gamma) * (t/2)^gamma
#' # it also holds for a complex variable:
#' t <- 1 + 2i
#' Bessel::BesselJ(t, gamma)
#' Bessel(m=15, t^2/4, gamma) * (t/2)^gamma
Bessel <- function(m, x, gamma){
  stopifnot(is.numeric(gamma) || is.complex(gamma), Re(gamma) > -1)
  if(is.matrix(x)){
    p <- nrow(x)
  }else{
    p <- length(x)
  }
  hypergeomPFQ(m, NULL, gamma+(p+1)/2, -x) / mvgamma(gamma+(p+1)/2, p)
}
