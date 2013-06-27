rownames = letters[1:3]
colnames = LETTERS[1:3]
mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
int.mat = matrix(c(rep(1L,5),rep(2L,4)),ncol=3,dimnames=list(rownames,colnames))
levels = c("AA","BB")
char.mat = matrix()
  
back.dir = tempdir()
ds.data.file = file.path(back.dir,"bigmat","ds")
ds.data.file = normalizePath(ds.data.file)
ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
fs.data.file = file.path(back.dir,"bigmat","fs")
fs = BigMatrixFactor(int.mat,fs.data.file,3,3,list(rownames,colnames),levels=levels)

test_creation <- function() {
  back.dir = tempdir()
  
  ds = BigMatrix(mat,tempfile())
  fs = BigMatrixFactor(int.mat,tempfile(),levels=levels)
  checkTrue(validObject(ds))
  checkTrue(validObject(fs))
  checkEquals(ds[,],mat)
  checkEquals(as(fs[,],"matrix"),int.mat)
  
  ds = BigMatrix(ds$bigmat,tempfile())
  fs = BigMatrixFactor(fs$bigmat,tempfile(),levels=levels)
  checkTrue(validObject(ds))
  checkTrue(validObject(fs))
  checkEquals(ds[,],mat)
  checkEquals(as(fs[,],"matrix"),int.mat)

  ds = BigMatrix(x=NULL,tempfile(),3,3,list(rownames,colnames))
  ds[,] = mat
  fs = BigMatrixFactor(x=NULL,tempfile(),3,3,list(rownames,colnames),levels=levels)
  fs[,] = int.mat
  checkTrue(validObject( ds ))
  checkEquals(ds[,],mat)
  checkTrue(validObject( fs ))
  checkEquals(as(fs[,],"matrix"),int.mat)
  checkException( BigMatrix(x=1:4, tempfile()), "x must be NULL, scalar numeric, matrix, or bigmatrix.", silent=TRUE)
  checkEquals( BigMatrix(x=12, ncol=2, nrow=2, dimnames=list(LETTERS[1:2], LETTERS[1:2]), tempfile())[2, 2], 12, "BigMatrix creation with a scalar init value.")
  na_bm = BigMatrix(backingfile=tempfile(), nrow=3, ncol=3, dimnames=NULL)
  na_bmf = BigMatrixFactor(backingfile=tempfile(), nrow=3, ncol=3, dimnames=NULL, levels=LETTERS[1:3])
  checkTrue( all(is.na(na_bm[, ])), "Default init is to all NA for BigMatrix")
  checkTrue( all(is.na(na_bmf[, ])), "Default init is to all NA for BigMatrixFactor")

  checkException( BigMatrixFactor(7, backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4]), "Scalar init value for BigMatrixFactor should not be > length(levels)", silent=TRUE)
  checkException( BigMatrixFactor(0, backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4]), "Scalar init value for BigMatrixFactor should not be > length(levels)", silent=TRUE)
  checkException( BigMatrixFactor("X", backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4]), "Scalar character init value for BigMatrixFactor should be in levels", silent=TRUE)
  checkTrue( all(as.integer(BigMatrixFactor(3L, backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4])[, ]) == 3L), "Scalar init value for BigMatrixFactor should be in 1:length(levels)")
  checkTrue( all(as.integer(BigMatrixFactor(3, backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4])[, ]) == 3L), "Scalar init value for BigMatrixFactor should be in 1:length(levels), double OK")
  checkTrue( all(BigMatrixFactor("B", backingfile=tempfile(), ncol=3, nrow=4, levels=LETTERS[1:4])[, ] == "B"), "Scalar character init value in levels is OK.")
}

test_coercion <- function() {
  checkEquals( ds[,], as(ds,"matrix") )
  checkEquals( ds[,], as.matrix(ds) )
  checkEquals( as(fs[,],"matrix"), as(fs,"matrix") )
  checkEquals( fs[,], as(fs,"factor") )
  checkEquals( fs[,], as.matrix(fs) )
}

test_subset <- function() {
  checkIdentical( ds[1,3], mat[1,3] )
  checkIdentical( ds[1,], mat[1,] )
  checkIdentical( ds[,2], mat[,2] )
  checkIdentical( ds[,], mat )
  checkIdentical( ds["a",], mat["a",] )
  checkIdentical( ds[,"B"], mat[,"B"] )
  ds["b","B"] = 3
  mat["b","B"] = 3
  checkIdentical( ds[,], mat,"After re-setting some values" )
  char.mat = matrix(c(rep(1L,5),rep(2L,4)),ncol=3,dimnames=dimnames(fs))
  char.mat = structure(char.mat, levels=levels(fs), class="factor")
  checkIdentical( fs[,], char.mat )
  checkIdentical( fs[,1], char.mat[,1] )
  checkIdentical( fs[2,], char.mat[2,] )
}

test_write <- function() {
  ds[1,1] = 5
  checkIdentical(ds[1,1],5,"Writing to a BigMatrix")
  ds[2,] = 5
  checkIdentical(ds[2,],c(A=5,B=5,C=5),"Writing row to a BigMatrix")
  ds[,2] = 7
  checkIdentical(ds[,2],c(a=7,b=7,c=7),"Writing col to a BigMatrix")
  ds[,] = mat
  checkIdentical(ds[,],mat,"Writing whole matrix to a BigMatrix")
  Sys.chmod(ds$datapath,"0444")
  ds$attach(force=TRUE)
  checkException( ds[1,1] <- 5, silent=TRUE, "Writing to a BigMatrix with a non-writeable data file")
  Sys.chmod(ds$datapath,"0644")
  ds$attach(force=TRUE)
  fs[,] = int.mat
  fs[,1] = 1:3
  returned.factor = factor(structure(c("AA","BB",NA),names=letters[1:3]),levels=levels)
  checkEquals(fs[,1], returned.factor, "Setting BigMatrixFactor with integers")
  fs[,1] = c("BB",NA,"AA")
  returned.factor = factor(structure(c("BB",NA,"AA"),names=letters[1:3]),levels=levels)
  checkIdentical(fs[,1], returned.factor, "Setting BigMatrixFactor with characters")
}

test_describing <- function() {
  mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
  ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
  # Getting
  checkEquals(nrow(ds),nrow(mat))
  checkEquals(ncol(ds),ncol(mat))
  checkEquals(dim(ds),dim(mat))
  checkEquals(length(ds),length(mat))
  checkIdentical(dimnames(ds),dimnames(mat))
  # Setting
  new.dimnames = list(letters[4:6], LETTERS[4:6])
  dimnames(ds) = new.dimnames
  checkIdentical(dimnames(ds), new.dimnames)
  
  ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
  colnames(ds) = new.dimnames[[2]]
  checkIdentical(colnames(ds), new.dimnames[[2]])

  ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
  rownames(ds) = new.dimnames[[1]]
  checkIdentical(rownames(ds), new.dimnames[[1]])
}

test_reattach <- function() {
  mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
  data.file = tempfile()
  old.ds = BigMatrix(mat,data.file,3,3,list(rownames,colnames))
  object.file = paste(data.file,".rds",sep="")
  saveRDS(old.ds,file=object.file)
  new.ds = readRDS(object.file)
  checkEquals( new.ds[,], mat )

  if (Sys.info()[["sysname"]] != "Windows") {
    Sys.chmod(old.ds$descpath,"0000")
    checkException(old.ds$attach(force=TRUE),silent=TRUE)
    Sys.chmod(old.ds$descpath,"0644")
  }
  
  file.rename(data.file,file.path(tempdir(),"shoe"))
  Sys.chmod(old.ds$descpath,"0600")
  if (Sys.info()[["sysname"]] != "Windows") {
    checkException(old.ds$attach(force=TRUE),silent=TRUE,"Missing backing file")
  }
  file.rename(old.ds$descpath,file.path(tempdir(),"shoe2"))
  checkException(old.ds$attach(force=TRUE),silent=TRUE,"Missing descriptor file")
}

test_levels <- function() {
  checkIdentical( levels(fs), c("AA","BB") )
  checkIdentical( fs$levels, c("AA","BB") )
  checkException( fs$levels <- c("AB","QQ"), silent=TRUE )
  checkException( levels(fs) <- c("AB","QQ"), silent=TRUE )
}

test_nlevels <- function() {
  checkEquals( nlevels(fs), 2L )
}

test_paths <- function() {
  desc.file = paste(ds.data.file,"desc.rds",sep=".")
  checkEquals(normalizePath(ds$datapath),ds.data.file)
  checkEquals(normalizePath(ds$descpath),desc.file)
  new.data.file = tempfile()
  new.ds = BigMatrix(mat,new.data.file,3,3,list(rownames,colnames))
  even.newer.data.file = tempfile()
  checkException({new.ds$datapath = even.newer.data.file},silent=TRUE,"backingfile must exist before datapath is replaced.")
  file.copy(ds.data.file,even.newer.data.file)
  new.ds$datapath = even.newer.data.file
  even.newer.data.file = normalizePath(even.newer.data.file)
  checkEquals(normalizePath(new.ds$datapath),even.newer.data.file)
  checkEquals(new.ds[,], ds[,])
}

test_apply <- function() {
  checkEquals( apply( ds, 1, mean ), apply( mat, 1, mean ), "apply on BigMatrix just like apply on base matrix" )
}

