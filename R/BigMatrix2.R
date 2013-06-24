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
##' @importClassesFrom Biobase AnnotatedDataFrame AssayData
##' @importFrom Biobase annotatedDataFrameFrom assayDataElementNames assayDataElement
##' @importFrom BiocGenerics updateObject
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

BigMatrix2Generator <- setRefClass("BigMatrix2",
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
                           backingfile="character",
                           .rownames="character", 
                           .colnames="character",
                           description="list",
                           .bm="big.matrix"
                           ),
                         methods=list(
                           rownames = function(value) {
                             if (missing(value)) {
                               return(.self$.rownames)
                             } else 
                               if (length(value) != .self$nrow) { stop("Length of rownames must match the number of rows.") }
                             .self$.rownames = value
                           }, 
                           colnames = function(value) {
                             if (missing(value)) {
                               return(.self$.colnames)
                             } else
                               if (length(value) != .self$nncol) { stop("Length of colnames must match the number of columns.") }
                             .self$.colnames = value
                           }, 
                           dimnames = function(value) {
                             if (missing(value)) {
                               return(list(.self$.rownames, .self$.colnames))
                             } else {
                               if (length(value) != 2) { stop("dimnames must be set with a list of length 2.") }
                               .self$.rownames(value[[1]])
                               .self$.colnames(value[[2]])
                             }
                           },
                           dim  = function() {
                             return( c( .self$description$nrow, .self$description$ncol ) )
                           },
                           nrow = function() {
                             return(.self$description$nrow)
                           },
                           ncol = function() {
                             return(.self$description$ncol)
                           },
                           length = function() {
                             return( .self$description$nrow * .self$description$ncol )
                           },
                           attach=function(force=FALSE) {
                             if (force == FALSE && ! is.nil(.self$.bm@address)) {
                               message("Already attached to on-disk data. To re-attach, use force=TRUE.\n")
                             } else {
                               message("Attaching to on-disk data:", .self$backingfile, "...\n")
                               if ( ! file.exists(backingfile) ) {
                                 stop("Backing file ",backingfile," does not exist.")
                               }
                               if ( file.access(.self$backingfile,4) != 0 ) {
                                 stop("Can not attach to descriptor file without read permissions.")
                               }
                               tryCatch({
                                 desc = new('big.matrix.descriptor',  description=.self$description)
                                 .self$.bm = attach.big.matrix(desc,path=dirname(.self$backingfile))
                               },
                                        error = function(e) { stop("Failed to attach big.matrix on disk component.\n") } )
                             }
                             if ("datapath" %in% ls(.self)) { warning("Attaching an older type of BigMatrix. Use updateObject method to update.") }
                           },
                           getValues=function(i,j,drop=TRUE, withDimnames=TRUE) {
                             object = .self$bigmat
                             if (!missing(i) && is.character(i)) { i = match(i,.self$.rownames) }
                             if (!missing(j) && is.character(j)) { j = match(j,.self$.colnames) }
                             if (missing(i)) {
                               if (missing(j)) {
                                 x = object[,,drop=drop]
                               } else {
                                 x = object[,j,drop=drop]
                                 }
                             } else {
                               if (missing(j)) {
                                 x = object[i,,drop=drop]
                               } else {
                                 x = object[i,j,drop=drop]
                               }
                             }
                             if (withDimnames == TRUE && base::length(x) > 1) {
                               if (is.matrix(x)) {
                                 if (missing(i)) { i = seq.int(1, .self$nrow())}
                                 if (missing(j)) { j = seq.int(1, .self$ncol())}
                                 dimnames(x) = list(.self$.rownames[i], .self$.colnames[j])
                               } else {
                                 if (missing(i)) { names(x) = .self$.rownames }
                                 if (missing(j)) { names(x) = .self$.colnames }
                               }
                             }
                             return(x)
                           },
                           setValues=function(i,j,value) {
                             object = .self$bigmat
                             bigmemory:::checkReadOnly(object)
                             if (!missing(i) && is.character(i)) { i = match(i,.self$.rownames) }
                             if (!missing(j) && is.character(j)) { j = match(j,.self$.colnames) }
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
                           save=function(rdsfile) {
                             saveRDS( .self, rdsfile )
                           },
                           show=function() {
                             message( class(.self), "\nbackingfile :", .self$backingfile, "\ndim: ", paste(.self$dim(), collapse=", "), "\n")
                             if (is.nil(.self$.bm@address)) {
                               message("Object is not currently attached to on-disk data.\n")
                             } else {
                               message("Object is currently attached to on-disk data.\n")
                             }
                           }
                           )
                         )

setValidity("BigMatrix2", function(object) {
  if (!file.exists(object$backingfile)) {
    return("Backingfile does not exist")
  }
  if (file.access(object$backingfile,4) < 0) {
    return("Backingfile is not readable")
  }
  return(TRUE)
})

## Maintain Generic Function and Method illusion for certain matrix API functions
## No logic other than passing to the right R5 method

##' @exportMethod '['
setMethod('[', signature(x = "BigMatrix2"),
          function(x,i,j,..., drop=TRUE) {
            return(x$getValues(i,j,..., drop=drop))
          })

##' @exportMethod '[<-'
setMethod('[<-', signature(x = "BigMatrix2",i="ANY",j="ANY",value="ANY"),
          function(x,i,j,value) {
            x$setValues(i,j,value)
            return(x)
          })

##' @exportMethod dimnames
setMethod('dimnames', signature(x="BigMatrix2"),
          function(x) {
            return(x$dimnames())
          })

##' @exportMethod 'dimnames<-'
setMethod('dimnames<-', signature(x="BigMatrix2",value="ANY"),
          function(x,value) {
            x$dimnames(value)
            return(x)
          })

##' @exportMethod nrow
setMethod('nrow', signature(x="BigMatrix2"),
          function(x) {
            return(x$nrow())
          })

##' @exportMethod ncol
setMethod('ncol', signature(x="BigMatrix2"),
          function(x) {
            return(x$ncol())
          })

##' @exportMethod dim
setMethod('dim', signature(x="BigMatrix2"),
          function(x) {
            return(x$dim())
          })

##' @exportMethod length
setMethod('length', signature(x="BigMatrix2"),
          function(x) {
            return(x$length())
          })

##' @exportMethod as.matrix
setMethod("as.matrix",signature(x="BigMatrix2"), function(x) { return(x[,]) })
setAs("BigMatrix2","matrix", function(from) { return(from[,]) })

##' @exportMethod apply
setMethod("apply",signature(X="BigMatrix2"), function(X, MARGIN, FUN, ...) { apply(X$bigmat, MARGIN, FUN, ...) })

##' Create a new BigMatrix2-derived class
##'
##' Create a new BigMatrix2-derived class
##' @param x NULL, matrix, or big.matrix. Optional data or big.matrix for new BigMatrix
##' @param class character, class name.  BigMatrix or BigMatrixFactor currently.
##' @param backingfile character, full path to the file that will contain the data matrix
##' @param nrow integer, number of rows in the matrix we are about to create
##' @param ncol integer, number of columns in the matrix we are about to create
##' @param dimnames list, list(rownames,colnames), as for a typical matrix
##' @param type character, can be double, integer, or char
##' @param ... other args to pass to "new" method of generator object, like levels
##' @return BigMatrix2
##' @keywords internal
##' @rdname initBigMatrix2
.initBigMatrix2 = function(x=NULL, class=c("BigMatrix2"), backingfile, nrow, ncol, dimnames=NULL, type="double", ...) {
  class = match.arg(class)
  backingpath = dirname(backingfile)
  dir.create(backingpath,showWarnings=FALSE,recursive=TRUE)
  backingpath = normalizePath(backingpath)
  descriptorfile = paste(backingfile,".desc",sep="")
  if (is.matrix(x)) {
    if (is.null(dimnames)) {
      dimnames = dimnames(x)
      dimnames(x) = NULL
    }
    new.matrix = as.big.matrix(x, backingpath=backingpath, descriptorfile=basename(descriptorfile), backingfile=basename(backingfile))
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
      backingfile=basename(backingfile),
      descriptorfile=descriptorfile,
      backingpath=backingpath)
  } else {
    stop("Argument x must be a scalar numeric, matrix, or big.matrix.\n")
  }
  unlink(file.path(backingpath,descriptorfile))  # Delete bigmemory version of desc file until they stop using dput/dget
  description = describe(new.matrix)@description # description method is not exported and it just does this anyway
  bm = getRefClass(class)$new(.bm=new.matrix, description=description, backingfile=backingfile, .rownames=dimnames[[1]], .colnames=dimnames[[2]], ...)
  if (!validObject(bm)) {
    stop("Failed to create a valid ", class, "!\n")
  }
  return( bm )
}

##' Create a new BigMatrix2
##'
##' Create a new BigMatrix2
##' @param x scalar numeric, NULL, matrix, or big.matrix. Optional data or big.matrix for new BigMatrix2. A scalar numeric can be used to initalize the whole matrix.
##' @param backingfile character, full path to the file that will contain the data matrix
##' @param nrow integer, number of rows in the matrix we are about to create
##' @param ncol integer, number of columns in the matrix we are about to create
##' @param dimnames list, list(rownames,colnames), as for a typical matrix
##' @param type character type of big.matrix (double, integer, char)
##' @return BigMatrix2
##' @examples
##' dnames = dimnames=list(letters[1:3],LETTERS[1:3])
##' x = matrix(1:9,ncol=3,dimnames=dnames)
##' ds = BigMatrix2(x,tempfile())
##' ds = BigMatrix2(backingfile=tempfile(),nrow=3,ncol=3,dimnames=dnames)
##' @export
BigMatrix2 <- function(x=NULL,backingfile,nrow,ncol,dimnames=NULL,type="double") {
  bm = .initBigMatrix2(x=x, class="BigMatrix2",backingfile=backingfile, nrow=nrow, ncol=ncol, dimnames=dimnames, type=type)
  return( bm )
}

##' Update previous BigMatrix type to new BigMatrix
##'
##' BigMatrix has changed some of its internal storage to eliminate the descriptor file and to keep the dimnames on the R side. This function will take a
##' Version <= 1.3 type and update it to the Version >=1.4 type.
##' @param object BigMatrix
##' @export 
##' @return BigMatrix
setMethod("updateObject", signature=signature(object="BigMatrix"), function(object) {
  tryCatch(
    { desc = readRDS(object$descpath) }, 
    error = function(e) { stop(sprintf("Failed to read descriptor file when updating BigMatrix.%s\n", e))} )
  desc.list = desc@description
  dimnames = list(desc.list$rowNames, desc.list$colNames)
  desc.list$rowNames = desc.list$colNames = NULL
  new.desc = new('big.matrix.descriptor',  description=desc.list)
  bm = attach.resource(new.desc, path=dirname(object$datapath))
  if (class(object) == "BigMatrixFactor") {
    bigmat = BigMatrix2(x=bm, backingfile=object$datapath, dimnames=dimnames, levels=object$levels)
  } else {
    bigmat = BigMatrix2(x=bm, backingfile=object$datapath, dimnames=dimnames)
  }
  return(bigmat)
})

# benchmark(x = theta$getValues(,  1),  y = theta2$getValues(,  1),  z = theta2$getValues(,  1,  withDimnames=FALSE),  replications=5)
