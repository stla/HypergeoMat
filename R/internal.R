#' @importFrom arrangements npartitions
#' @importFrom utils head tail
NULL

.P <- function(m, n){
  sum(sapply(1L:m, function(i){
    sum(sapply(1L:min(m,n,i), function(k){
      npartitions(i, k)
    }))
  }))
}

.T <- function(alpha, a, b, kappa, i){
  if(length(kappa) == 0L || kappa[1L] == 0L) return(1)
  c <- kappa[i] - 1L - (i-1L)/alpha
  d <- kappa[i]*alpha - i
  s <- seq_len(kappa[i]-1L)
  e <- d - s*alpha + vapply(s, function(i) sum(kappa >= i), integer(1L))
  g <- e + 1
  s <- seq_len(i-1L)
  f <- kappa[s]*alpha - s - d
  h <- f + alpha
  l <- h*f
  prod1 <- prod(a+c)
  prod2 <- prod(b+c)
  prod3 <- prod((g-alpha)*e / (g*(e+alpha)))
  prod4 <- prod((l-f) / (l+h))
  out <- prod1/prod2 * prod3 * prod4
  return(ifelse(is.nan(out) || is.infinite(out), 0, out))
} # prod3 prod4 NaN pour alpha = Inf ; T -> 0; out*alpha -> prod1/prod2 / 2

.betaratio <- function(kappa, mu, k, alpha){
  t <- k - alpha*mu[k]
  s <- seq_len(k)
  u <- t + 1 - s + alpha*kappa[s]
  s <- seq_len(k-1L)
  v <- t - s + alpha*mu[s]
  s <- seq_len(mu[k]-1L)
  w <- vapply(s, function(i) sum(mu >= i), integer(1L)) - t - alpha*s
  alpha * prod(u/(u+alpha-1)) * prod((v+alpha)/v) * prod((w+alpha)/w)
} # NaN pour alpha = Inf ; tend vers Inf
  # 1er prod -> 0.5 (dépend de kappa et de mu[k]), 2ème prod -> 1; 3ème prod -> 4 (kappa=(5,2) mu = (4,1) et k=1)
  # 3ème prod semble être mu[k]

HypergeoI <- function(m, alpha, a, b, n, x){
  summation <- function(i, z, j){
    for(kappai in seq_len(if(i==1L) j else min(kappa[i-1L],j))){
      kappa[i] <<- kappai
      t <- .T(alpha, a, b, kappa, i)
      z <- z * x * (n-i+1L+alpha*(kappai-1L)) * t
      s <<- s + z
      if(j > kappai && i < n){
        summation(i+1L, z, j-kappai)
      }
    }
    kappa[i] <<- 0L
  }
  s <- 1
  kappa <- integer(0L)
  summation(1L, 1, m)
  s
}

isPositiveInteger <- function(m){
  is.vector(m) && is.numeric(m) && length(m) == 1L && floor(m) == m
}

isSymmetricPositive <- function(M){
  isSymmetric(M) && all(eigen(M, symmetric = TRUE, only.values = TRUE)$values >= 0)
}

isNotNegativeInteger <- function(z){
  Im(z) != 0 || Re(z)>0 || Re(z) != trunc(Re(z))
}

isNumericOrComplex <- function(x){
  is.vector(x) && is.atomic(x) && (is.numeric(x) || is.complex(x))
}
