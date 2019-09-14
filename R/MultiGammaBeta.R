#' Multivariate Gamma function (of complex variable)
#'
#' @description The multivariate Gamma function (\code{mvgamma}) and
#' its logarithm (\code{lmvgamma}).
#'
#' @name mvgamma
#'
#' @param x a real or a complex number with \code{Re(x)>0}
#' @param p a positive integer, the dimension
#'
#' @return A real or a complex number.
#' @importFrom gsl lngamma_complex
#'
#' @examples x <- 5
#' mvgamma(x, p = 2)
#' sqrt(pi)*gamma(x)*gamma(x-1/2)
NULL

.lmvgamma <- function(x, p){
  C <- p*(p-1)/4*log(pi)
  if(is.numeric(x)){
    S <- sum(lgamma(x + (1L - seq_len(p))/2))
  }else{
    S <- sum(lngamma_complex(Re(x) + (1L - seq_len(p))/2), rep(Im(x), p))
  }
  C + S
}


#' @rdname mvgamma
#' @export
lmvgamma <- function(x, p){
  stopifnot(
    isPositiveInteger(p),
    is.numeric(x) || is.complex(x),
    length(x) == 1L,
    Re(x) > 0
  )
  .lmvgamma(x,p)
}

#' @rdname mvgamma
#' @export
mvgamma <- function(x, p){
  exp(lmvgamma(x, p))
}


#' Multivariate Beta function (of complex variable)
#'
#' @description The multivariate Beta function (\code{mvbeta}) and
#' its logarithm (\code{lmvbeta}).
#'
#' @name mvbeta
#'
#' @param a,b real or complex numbers with \code{Re(a)>0}, \code{Re(b)>0}
#' @param p a positive integer, the dimension
#'
#' @return A real or a complex number.
#'
#' @examples a <- 5; b <- 4; p <- 3
#' mvbeta(a, b, p)
#' mvgamma(a, p) * mvgamma(b, p) / mvgamma(a+b, p)
NULL

.lmvbeta <- function(a, b, p){
  .lmvgamma(a, p) + .lmvgamma(b, p) - .lmvgamma(a+b, p)
}

#' @rdname mvbeta
#' @export
lmvbeta <- function(a, b, p){
  stopifnot(
    isPositiveInteger(p),
    is.numeric(a) || is.complex(a),
    length(a) == 1L,
    Re(a) > 0,
    is.numeric(b) || is.complex(b),
    length(b) == 1L,
    Re(b) > 0
  )
  .lmvbeta(a, b, p)
}

#' @rdname mvbeta
#' @export
mvbeta <- function(a, b, p){
  exp(lmvbeta(a, b, p))
}
