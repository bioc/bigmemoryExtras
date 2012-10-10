.onLoad <- function(libname,pkgname) {
  options(bigmemory.typecast.warning=FALSE)
#  options(bigmemory.allow.dimnames=TRUE)
}

.onUnLoad <- function(libpath) {
  options(bigmemory.typecast.warning=NULL)
#  options(bigmemory.allow.dimnames=FALSE)
}
