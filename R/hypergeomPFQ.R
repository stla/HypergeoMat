#' Hypergeometric function of a matrix argument
#'
#' @description Evaluates a truncated hypergeometric function of a matrix
#' argument.
#'
#' @param m truncation weight, a positive integer
#' @param a the "upper" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param b the "lower" parameters, a numeric or complex vector,
#' possibly empty (or \code{NULL})
#' @param x either a real or complex square matrix,
#' or a numeric or complex vector, the eigenvalues of the matrix
#' @param alpha the alpha parameter, a positive number
#'
#' @return A real or a complex number.
#' @export
#'
#' @note The hypergeometric function of a matrix argument is usually defined
#' for a symmetric real matrix or a Hermitian complex matrix.
#'
#' @references Koev, P. and Edelman, A. (2006).
#' \emph{The Efficient Evaluation of the Hypergeometric Function of a Matrix Argument}.
#' Mathematics of Computation, 75, 833-846, 2006.
#'
#' @examples # a scalar x example, the Gauss hypergeometric function
#' hypergeomPFQ(m = 20, a = c(1,2), b = c(3), x = 0.5)
#' gsl::hyperg_2F1(1, 2, 3, 0.5)
#' # 0F0 is the exponential of the trace
#' X <- toeplitz(c(3,2,1))/10
#' hypergeomPFQ(m = 10, a = NULL, b = NULL, x = X)
#' exp(sum(diag(X)))
#' # 1F0 is det(I-X)^(-a)
#' X <- toeplitz(c(3,2,1))/100
#' hypergeomPFQ(m = 15, a = 3, b = NULL, x = X)
#' det(diag(3)-X)^(-3)
#' # Herz's relation for 1F1
#' hypergeomPFQ(m=15, a = 2, b = 3, x = X)
#' exp(sum(diag(X))) * hypergeomPFQ(m=15, a = 3-2, b = 3, x = -X)
#' # Herz's relation for 2F1
#' hypergeomPFQ(15, a = c(1,2), b = 3, x = X)
#' det(diag(3)-X)^(-2) *
#'   hypergeomPFQ(15, a = c(3-1,2), b = 3, -X%*%solve(diag(3)-X))
hypergeomPFQ <- function(m, a, b, x, alpha = 2){
  stopifnot(
    isPositiveInteger(m),
    is.null(a) || is.atomic(a),
    is.null(b) || is.atomic(b),
    is.atomic(alpha),
    length(alpha) == 1L,
    is.numeric(alpha),
    alpha > 0
  )
  if(is.matrix(x)){
    # stopifnot(isSymmetric(x))
    # x <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    x <- eigen(x, only.values = TRUE)$values
  }else{
    stopifnot(is.atomic(x), is.numeric(x) || is.complex(x))
  }
  if(all(x == x[1L])){
    return(HypergeoI(m, alpha, a, b, length(x), x[1L]))
  }
  #
  summation <- function(i, z, j){
    jack <- function(k, beta, c, t, mu, Nmu){
      for(i in max(k,1L) : sum(mu > 0L)){
        if(length(mu) == i || mu[i] > mu[i+1L]){
          d <- Nmu
          gamma <- beta * .betaratio(kappa, mu, i, alpha)
          mu[i] <- mu[i]-1L
          Nmu <- .Nkappa(mu) - 1L
          if(mu[i] > 0L){
            jack(i, gamma, c+1L, t, mu, Nmu)
          }else{
            if(Nkappa > 1L){
              J[Nkappa,t] <<- J[Nkappa,t] + gamma*x[t]^(c+1L) *
                ifelse(mu[1]>0L, J[Nmu,t-1L], 1) #any(mu>0L) <=> mu[1]>0L ?
            }
          }
          mu[i] <- mu[i] + 1L
          Nmu <- d
        }
      }
      if(k == 0L){
        if(Nkappa>1L) J[Nkappa,t] <<- J[Nkappa,t] + J[Nkappa,t-1L]
      }else{
        J[Nkappa,t] <<- J[Nkappa,t] + beta * x[t]^c * J[Nmu,t-1L]
      }
    }
    # end jack
    r <- if(i == 1L) j else min(kappa[i-1L],j)
    for(kappai in seq_len(r)){
      kappa[i] <<- kappai
      Nkappa <- .Nkappa(kappa) - 1L
      z <- z * .T(alpha, a, b, kappa, i)
      if(Nkappa > 1L && (length(kappa) == 1L || kappa[2L] == 0L)){
        J[Nkappa,1L] <<- x[1L]*(1+alpha*(kappa[1L]-1L)) * J[Nkappa-1L,1L]
      }
      for(t in 2L:n) jack(0L, 1, 0L, t, kappa, Nkappa)
      s <<- s + z * J[Nkappa,n]
      if(j > kappai && i < n){
        summation(i+1L, z, j-kappai)
      }
    }
    kappa[i] <<- 0L
  } # end summation
  #
  #
  n <- length(x)
  Pmn <- .P(m,n)
  s <- 1
  J <- matrix(0, Pmn, n)
  J[1L,] <- cumsum(x)
  kappa <- integer(0L)
  #
  D <- rep(NA_integer_, Pmn)
  Last <- t(as.integer(c(0,m,m)))
  fin <- 0L
  for(. in 1L:n){
    NewLast <- matrix(NA_integer_, nrow = 0L, ncol = 3L)
    for(i in 1L:nrow(Last)){
      manque <- Last[i,2L]
      l <- min(manque, Last[i,3L])
      if(l > 0L){
        D[Last[i,1L]+1L] <- fin + 1L
        s <- 1L:l
        NewLast <- rbind(NewLast, cbind(fin+s, manque-s, s))
        fin <- fin + l
      }
    }
    Last <- NewLast
  }
  .Nkappa <- function(kappa){
    kappa <- kappa[kappa>0L]
    if(length(kappa) == 0L) return(1L)
    D[.Nkappa(head(kappa,-1L))] + tail(kappa,1L)
  }
  summation(1L, 1, m)
  s
}
