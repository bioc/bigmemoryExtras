##' An extension of the bigmemory package with added safety, convenience, and a factor class.
##'
##' This package defines a "BigMatrix" ReferenceClass which adds
##' safety and convenience features to the filebacked.big.matrix class from the
##' bigmemory package. BigMatrix protects against segfaults by monitoring, and 
##' gracefully restoring, the connection to on-disk data. It also protects 
##' against accidental data modification with a filesystem-based permissions 
##' system. We provide functionality for using BigMatrix-derived classes 
##' as assayData matrices within the Biobase package's eSet family of classes. 
##' BigMatrix provides some optimizations related to attaching to,
##' and indexing into, file-backed matrices with dimnames. Additionally, the package 
##' provides a "BigMatrixFactor" class, a file-backed matrix with factor 
##' properties.
##'
##' BigMatrix stores the
##' filesystem path to its on-disk data. When a big.matrix object is saved, e.g. with
##' "save" the external pointer to the on-disk data becomes "nil". Any access
##' to this nil pointer causes a segfault. A BigMatrix object will
##' gracefully re-attach itself to its on-disk components any time use of the
##' "nil" pointer is attempted.
##'
##' @docType package
##' @name bigmemoryExtras-package
##' @aliases bigmemoryExtras-package
##' @seealso BigMatrix-class BigMatrixFactor-class BigMatrix BigMatrixFactor filebacked.big.matrix ReferenceClasses
##' @import methods
##' @import bigmemory
##' @importFrom biganalytics apply
NULL

### Modifications to big.matrix object from bigmemory to help with saving to and restoring from disk and to prevent usage of a nil address

# What to do with description?  Currently really slow to attach as reading dimnames from desc file with dget is 10X the time of reading same info with load
# Option 1: Dimnames and description as fields, use desc file only for checking read permissions.
#    Would need to overwrite getters and setters: "[", dimnames, object creation would need to put these in their slots too
# Option 2: Replace description file with a RData file with same content.  Use finalize method to flush dimnames to desc file. Dimnames as field, or let big.matrix keep them?
#    Dimnames stay with data file +/- sync errors.  Data files useful independent of object. Maybe that is bad as col/row ordering for one element of a GenoSet could change.
#    Would be simpler overall
#    Would need save, finalize methods
#    Who does a better job using dimnames, R or big.matrix? Speed of indexing, etc.
#    I'll have to read the desc myself, so should I keep it in a field and make use of it?  Remove dimnames from desc for attach to avoid giving them to big.matrix?
#    Ooh, could use readRDS and saveRDS, like load/save, but for one object.  No need to deal with getting names of stuff from load.
#    Maybe want third file with rds version of desc so it could still be a useful big.matrix object?
#    Could save dimnames in finalize, but useful to have many objects pointing to the same big.matrix (SNOW, etc.).  Don't want contention. Don't have locking?

# Initialize can't require args or have side-effects, so can not create on-disk elements with initialize. Have to do "BigMatrix" function (and/or method?) for users to initialize

########################
###  Class BigMatrix ###
########################

BigMatrixGenerator <- setRefClass("BigMatrix",
                         fields=list(
                           bigmat=function(value) {
                             if (missing(value)) {
                               if (is.nil(.self$.bm@address)) {
                                 .self$attach()
                               }
                               return(.self$.bm)
                             } else {
                               if (is.big.matrix(value)) {
                                 .self$.bm = value
                               } else {
                                 stop("Replacement value for bigmat must be a big.matrix\n")
                               }
                             }
                           },
                           datapath=function(value) {
                             desc = describe(.self$bigmat)
                             if (missing(value)) {
                               backingfile = file.path(dirname(.self$descpath),desc@description$filename)
                               return(backingfile)
                             } else {
                               if (!file.exists(value)) {
                                 stop("Replacement datapath file must exist before setting 'datapath'. We will need to use the file to reattach after replacing datapath.\n")
                               }
                               desc@description$filename = basename(value)
                               saveRDS(desc,file=.self$descpath)
                               .self$attach(force=TRUE)
                               validObject(.self)
                             }
                           },
                           descpath="character",
                           .bm="big.matrix"
                           ),
                         methods=list(
                           colnames = function(value) {
                             if (missing(value)) {
                               return(base::colnames(.self$bigmat))
                             }
                             base::colnames(.self$bigmat) = value
                           },
                           rownames = function(value) {
                             if (missing(value)) {
                               return(base::rownames(.self$bigmat))
                             }
                             base::rownames(.self$bigmat) = value
                           },
                           dimnames = function(value) {
                             if (missing(value)) {
                               return(base::dimnames(.self$bigmat))
                             }
                             base::dimnames(.self$bigmat) = value
                           },
                           dim  = function() {
                             return(base::dim(.self$bigmat))
                           },
                           nrow = function() {
                             return(base::nrow(.self$bigmat))
                           },
                           ncol = function() {
                             return(base::ncol(.self$bigmat))
                           },
                           length = function() {
                             return(base::length(.self$bigmat))
                           },
                           attach=function(force=FALSE) {
                             if (force == FALSE && ! is.nil(.self$.bm@address)) {
                               message("Already attached to on-disk data. To re-attach, use force=TRUE.\n")
                             } else {
                               message("Attaching to on-disk data:", .self$descpath, "...\n")
                               if ( ! file.exists(.self$descpath) ) {
                                 stop("Descriptor file ",.self$descpath," does not exist.")
                               }
                               if ( file.access(.self$descpath,4) != 0 ) {
                                 stop("Can not attach to descriptor file without read permissions.")
                               }
                               tryCatch( { desc = readRDS(.self$descpath) },
                                        error = function(e) {  simpleError("Failed to attach big.matrix on disk component.\n") } )
                               backingfile = file.path(dirname(.self$descpath),desc@description$filename)
                               if ( ! file.exists(backingfile) ) {
                                 stop("Backing file ",backingfile," does not exist.")
                               }
                               if ( file.access( backingfile, 4 ) != 0 ) {
                                 stop("Can not attach to backing file without read permissions on the backing file.")
                               }
                               tryCatch({
                                 .self$.bm = attach.big.matrix(desc,path=dirname(.self$descpath))
                               },
                                        error = function(e) {  simpleError("Failed to attach big.matrix on disk component.\n") } )
                             }
                           },
                           getValues=function(i,j,drop=TRUE) {
                             object = .self$bigmat
                             if (!missing(i) && is.character(i)) { i = match(i,base::rownames(object)) }
                             if (!missing(j) && is.character(j)) { j = match(j,base::colnames(object)) }
                             if (missing(i)) {
                               if (missing(j)) {
                                 return(object[,,drop=drop])
                               } else {
                                 return(object[,j,drop=drop])
                               }
                             } else {
                               if (missing(j)) {
                                 return(object[i,,drop=drop])
                               } else {
                                 return(object[i,j,drop=drop])
                               }
                             }
                           },
                           setValues=function(i,j,value) {
                             object = .self$bigmat
                             bigmemory:::checkReadOnly(object)
                             if (!missing(i) && is.character(i)) { i = match(i,base::rownames(object)) }
                             if (!missing(j) && is.character(j)) { j = match(j,base::colnames(object)) }
                             if (missing(i)) {
                               if (missing(j)) {
                                 object[,] <- value
                               } else {
                                 object[,j] <- value
                                 }
                             } else {
                               if (missing(j)) {
                                 object[i,] <- value
                               } else {
                                 object[i,j] <- value
                               }
                             }
                           },
                           save=function() {
                             saveRDS( describe(.self$bigmat), file=.self$descpath )
                           },
                           show=function() {
                             message( class(.self), "\ndescpath:", .self$descpath, "\ndatapath:", .self$datapath, "\nnrow:", .self$nrow(), "\nncol:", .self$ncol(), "\n")
                             if (is.nil(.self$.bm@address)) {
                               message("Object is not currently attached to on-disk data.\n")
                             } else {
                               message("Object is currently attached to on-disk data.\n")
                             }
                           }
                           )
                         )

setValidity("BigMatrix", function(object) {
  if (!file.exists(object$descpath)) {
    return("Description file does not exist")
  } 
  if (!file.exists(object$datapath)) {
    return("Data file does not exist")
  }
  if (file.access(object$descpath,4) < 0) {
    return("Description file is not readable")
  }
  if (file.access(object$datapath,4) < 0) {
    return("Data file is not readable")
  }
  return(TRUE)
})

## Maintain Generic Function and Method illusion for certain matrix API functions
## No logic other than passing to the right R5 method

##' @exportMethod '['
setMethod('[', signature(x = "BigMatrix"),
          function(x,i,j,drop=TRUE) {
            return(x$getValues(i,j,drop=drop))
          })

##' @exportMethod '[<-'
setMethod('[<-', signature(x = "BigMatrix",i="ANY",j="ANY",value="ANY"),
          function(x,i,j,value) {
            x$setValues(i,j,value)
            return(x)
          })

##' @exportMethod dimnames
setMethod('dimnames', signature(x="BigMatrix"),
          function(x) {
            return(x$dimnames())
          })

##' @exportMethod 'dimnames<-'
setMethod('dimnames<-', signature(x="BigMatrix",value="ANY"),
          function(x,value) {
            x$dimnames(value)
            return(x)
          })

##' @exportMethod nrow
setMethod('nrow', signature(x="BigMatrix"),
          function(x) {
            return(x$nrow())
          })

##' @exportMethod ncol
setMethod('ncol', signature(x="BigMatrix"),
          function(x) {
            return(x$ncol())
          })

##' @exportMethod dim
setMethod('dim', signature(x="BigMatrix"),
          function(x) {
            return(x$dim())
          })

##' @exportMethod length
setMethod('length', signature(x="BigMatrix"),
          function(x) {
            return(x$length())
          })

##' @exportMethod as.matrix
setMethod("as.matrix",signature(x="BigMatrix"), function(x) { return(x[,]) })
setAs("BigMatrix","matrix", function(from) { return(from[,]) })

##' @exportMethod apply
setMethod("apply",signature(X="BigMatrix"), function(X, MARGIN, FUN, ...) { apply(X$bigmat, MARGIN, FUN, ...) })

##' Create a new BigMatrix-derived class
##'
##' Create a new BigMatrix-derived class
##' @param x NULL, matrix, or big.matrix. Optional data or big.matrix for new BigMatrix
##' @param class character, class name.  BigMatrix or BigMatrixFactor currently.
##' @param backingfile character, full path to the file that will contain the data matrix
##' @param nrow integer, number of rows in the matrix we are about to create
##' @param ncol integer, number of columns in the matrix we are about to create
##' @param dimnames list, list(rownames,colnames), as for a typical matrix
##' @param type character, can be double, integer, or char
##' @param ... other args to pass to "new" method of generator object, like levels
##' @return BigMatrix
##' @keywords internal
##' @rdname initBigMatrix
.initBigMatrix = function(x=NULL, class=c("BigMatrix","BigMatrixFactor"), backingfile, nrow, ncol, dimnames, type="double", ...) {
  class = match.arg(class)
  backingpath = dirname(backingfile)
  dir.create(backingpath,showWarnings=FALSE,recursive=TRUE)
  backingpath = normalizePath(backingpath)
  backingfile = basename(backingfile)
  descriptorfile = paste(backingfile,".desc",sep="")
  if (is.matrix(x)) {
    new.matrix = as.big.matrix(x, backingpath=backingpath, descriptorfile=descriptorfile, backingfile=backingfile)
  } else if (is.big.matrix(x)) {
    if( is.nil(x@address) ) {
      tryCatch( { x = attach.big.matrix(descpath) },
               error = function(e) { stop("Failed to attach big.matrix on disk component.\n") } )
    }
    new.matrix = x
  } else if (is.null(x) || (is.numeric(x) && length(x) == 1)) {
    new.matrix = filebacked.big.matrix(
      init=x, 
      nrow=nrow,ncol=ncol,
      type=type,
      backingfile=backingfile,
      descriptorfile=descriptorfile,
      backingpath=backingpath,
      dimnames=dimnames)
  } else {
    stop("Argument x must be a scalar numeric, matrix, or big.matrix.\n")
  }
  descpath = file.path(backingpath,paste(descriptorfile,".rds",sep=""))
  unlink(file.path(backingpath,descriptorfile))  # Delete bigmemory version of desc file until they stop using dput/dget
  bm = getRefClass(class)$new(.bm=new.matrix,descpath=descpath, ...)
  bm$save()
  if (!validObject(bm)) {
    stop("Failed to create a valid ", class, "!\n")
  }
  return( bm )
}

##' Create a new BigMatrix
##'
##' Create a new BigMatrix
##' @param x scalar numeric, NULL, matrix, or big.matrix. Optional data or big.matrix for new BigMatrix. A scalar numeric can be used to initalize the whole matrix. NULL gives the bigmemory feature of initializing to all zeros instantly.
##' @param backingfile character, full path to the file that will contain the data matrix
##' @param nrow integer, number of rows in the matrix we are about to create
##' @param ncol integer, number of columns in the matrix we are about to create
##' @param dimnames list, list(rownames,colnames), as for a typical matrix
##' @param type character type of big.matrix (double, integer, char)
##' @return BigMatrixFactor
##' @examples
##' dnames = dimnames=list(letters[1:3],LETTERS[1:3])
##' x = matrix(1:9,ncol=3,dimnames=dnames)
##' ds = BigMatrix(x,tempfile())
##' ds = BigMatrix(backingfile=tempfile(),nrow=3,ncol=3,dimnames=dnames)
##' @export
BigMatrix <- function(x=NA_real_,backingfile,nrow,ncol,dimnames=NULL,type="double") {
  bm = .initBigMatrix(x=x, class="BigMatrix",backingfile=backingfile, nrow=nrow, ncol=ncol, dimnames=dimnames, type=type)
  return( bm )
}
