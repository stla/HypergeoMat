#' @title Evaluation with Julia
#' @description Evaluate the hypergeometric function of a matrix argument with
#' Julia. This is highly faster.
#'
#' @return A function with the same arguments as \code{\link{hypergeomPFQ}}.
#'
#' @importFrom JuliaConnectoR juliaSetupOk juliaCall juliaImport
#' @export
hypergeomPFQ_julia <- function(){
  if(!juliaSetupOk()){
    stop("Julia setup is not OK.")
  }
  module <- system.file("julia", "HypergeoMat.jl", package = "HypergeoMat")
  . <- juliaCall("include", module)
  HypergeomPQ <- juliaImport(".HypergeomPQ", all = FALSE)
  function(m, a, b, x, alpha = 2){
    stopifnot(
      isPositiveInteger(m),
      is.null(a) || isNumericOrComplex(a),
      !anyNA(a),
      is.null(b) || isNumericOrComplex(b),
      !anyNA(b),
      isSquareMatrix(x) || isNumericOrComplex(x),
      length(x) != 0L,
      !anyNA(x),
      isNumber(alpha),
      alpha > 0
    )
    if(!is.matrix(x)){
      x <- unname(as.list(x))
    }
    HypergeomPQ$hypergeomPQ(
      as.integer(m), unname(as.list(a)), unname(as.list(b)), x, unname(alpha)
    )
  }
}
