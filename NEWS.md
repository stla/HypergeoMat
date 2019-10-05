# Version 2.0.0 (2019-10-05)

The package now uses `Rcpp` to evaluate the hypergeometric function when none 
of its arguments is a complex number.


# Version 1.0.1 (2019-09-26)

- Fixed `lmvgamma` for complex values

- Allows complex values `z` with `Re(z)<0` in `mvgamma`

- Added more unit tests


# Version 1.0.0 (2019-09-16)

First release.
