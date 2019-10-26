## Release summary

This is a major release. 
The package now uses `Rcpp` to evaluate the hypergeometric function in all cases, 
for real or complex values of the arguments, while the previous version used 
`Rcpp` only when the values of the arguments are real numbers.


## Test environments

* Ubuntu 16.04, R 3.4.4
* Windows 7 64-bit, R 3.6.1
* win-builder (devel and release)

## R CMD check results

Status: OK
