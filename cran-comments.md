This is a resubmission. The CRAN checks detected a C++ error. It was due to 
this line of code:

`int* ptr_data = &vout[0];`

which is incorrect when the vector `vout` is empty. So I've treated the empty 
case to fix this error.


## Testing environments

* Windows 10, R-4.1.2
* Ubuntu 18, R-3.6.3
* Ubuntu 20 via Github action
* win-builder devel
* mac-builder
* r-hub


## R CMD check results

OK
