
##############################
###  Class BigMatrixFactor ###
##############################

BigMatrixFactorGenerator <- setRefClass("BigMatrixFactor",
                             contains="BigMatrix",
                               fields=list(
                                 levels="character"
                                 ),
                               methods=list(
                                 getValues=function(i,j,drop=TRUE) {
                                   mat = callSuper(i,j,drop=drop)
                                   mat = structure(mat, levels=.self$levels, class="factor")
                                   return(mat)
                                 },
                                 setValues=function(i,j,value) {
                                   if (is.character(value)) {
                                     value = match(value,.self$levels)
                                   }
                                   callSuper(i,j,value)
                                 },
                                 nlevels=function() {
                                   return(base::length(.self$levels))
                                 },
                                 show=function() {
                                   callSuper()
                                   message("Levels:", paste(.self$levels, collapse=" "), "\n")
                                 }
                                )
                               )
BigMatrixFactorGenerator$lock("levels")

##' Create a new BigMatrixFactor
##'
##' Create a new BigMatrixFactor
##' @param x scalar numeric, NULL, matrix, or big.matrix. Optional data or big.matrix for new BigMatrix. A scalar numeric can be used to initalize the whole matrix. NULL gives the bigmemory feature of initializing to all zeros instantly.
##' @param backingfile character, full path to the file that will contain the data matrix
##' @param nrow integer, number of rows in the matrix we are about to create
##' @param ncol integer, number of columns in the matrix we are about to create
##' @param dimnames list, list(rownames,colnames), as for a typical matrix
##' @param levels character, as for a typical factor
##' @return BigMatrixFactor
##' @examples
##' dnames = dimnames=list(letters[1:3],LETTERS[1:3])
##' x = matrix( sample( 1:3, 9, replace=TRUE), ncol=3, dimnames=dnames)
##' ds = BigMatrixFactor(x,tempfile(),levels=c("AA","AB","BB"))
##' ds = BigMatrixFactor(backingfile=tempfile(),nrow=3,ncol=3,dimnames=dnames,levels=c("AA","AB","BB"))
##' @export
BigMatrixFactor <- function(x=NA_integer_,backingfile,nrow,ncol,dimnames,levels) {
  if ( length(levels) > 127 ) {
    type = "integer"
  } else {
    type = "char"
  }
  bm = .initBigMatrix(x=x,class="BigMatrixFactor",backingfile=backingfile, nrow=nrow, ncol=ncol, dimnames=dimnames, levels=levels, type=type)
  return( bm )
}

## Maintain Generic Function and Method illusion for certain matrix API functions
## No logic other than passing to the right R5 method
##' @exportMethod levels
setMethod("levels",signature=signature(x="BigMatrixFactor"),
          function(x) {
            return(x$levels)
          })

##' @exportMethod 'levels<-'
setMethod("levels<-",signature=signature(x="BigMatrixFactor"),
          function(x) {
            stop("Levels are read-only.")
          })

##' @exportMethod nlevels
setMethod("nlevels",signature=signature(x="BigMatrixFactor"),
          function(x) {
            return(x$nlevels())
          })

setAs("BigMatrixFactor","matrix", function(from) { return(from$bigmat[,]) })
setAs("BigMatrixFactor","factor", function(from) { return(from[,]) })
