library(bigmemoryExtras)
library(testthat)

rownames = letters[1:3]
colnames = LETTERS[1:3]
mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
int.mat = matrix(c(rep(1L,5),rep(2L,4)),ncol=3,dimnames=list(rownames,colnames))
levels = c("AA","BB")
char.mat = matrix()

back.dir = tempdir()
ds.data.file = file.path(back.dir,"bigmat","ds")
unlink(ds.data.file)
ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))

test_that("We can create BigMatrix objects", {
  back.dir = tempdir()
  back.file = tempfile()
  ds = BigMatrix(mat,back.file)
  expect_true(validObject(ds))
  expect_equal(ds[,],mat)

  ds = BigMatrix(ds$bigmat,back.file, dimnames=dimnames(mat))
  expect_true(validObject(ds))
  expect_equal(ds[,],mat)

  ds = BigMatrix(x=NULL,tempfile(),3,3,list(rownames,colnames))
  ds[,] = mat
  expect_true(validObject( ds ))
  expect_equal(ds[,],mat)
  expect_error( BigMatrix(x=1:4, tempfile()) )
  expect_equal( BigMatrix(x=12, ncol=2, nrow=2, dimnames=list(LETTERS[1:2], LETTERS[1:2]), tempfile())[2, 2], 12)
  na_bm = BigMatrix(backingfile=tempfile(), nrow=3, ncol=3, dimnames=NULL)
  expect_true( all(is.na(na_bm[, ])), "Default init is to all NA for BigMatrix")
})

test_that("We can coerce between types", {
  expect_equal( ds[,], as(ds,"matrix") )
  expect_equal( ds[,], as.matrix(ds) )
})

test_that("We can subset", {
  unlink(ds.data.file)
  ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
  expect_identical( ds[1,3], mat[1,3] )
  expect_identical( ds[1,], mat[1,] )
  expect_identical( ds[,2], mat[,2] )
  expect_identical( ds[,], mat )
  expect_identical( ds["a",], mat["a",] )
  expect_identical( ds[,"B"], mat[,"B"] )
  ds["b","B"] = 3
  mat["b","B"] = 3
  expect_identical( ds[,], mat,"After re-setting some values" )
  mat_nonames = mat
  dimnames(mat_nonames) = NULL
  ds_nonames = BigMatrix(mat_nonames,tempfile(),3,3)
  ds_nonames_return = ds_nonames[,];   dimnames(ds_nonames_return) = NULL  # Hack for bad handling of dimnames when they are NULL by bigmemory < 4.4.4
  expect_identical( ds_nonames_return, mat_nonames, "Get full matrix from BM w/o dimnames" )
})

test_that("We can write to objects", {
  ds[1,1] = 5
  expect_identical(ds[1,1],5,"Writing to a BigMatrix")
  ds[2,] = 5
  expect_identical(ds[2,],c(A=5,B=5,C=5),"Writing row to a BigMatrix")
  ds[,2] = 7
  expect_identical(ds[,2],c(a=7,b=7,c=7),"Writing col to a BigMatrix")
  ds[,] = mat
  expect_identical(ds[,],mat,"Writing whole matrix to a BigMatrix")
  Sys.chmod(ds$backingfile,"0444")
  ds$attach(force=TRUE)
  expect_error( ds[1,1] <- 5 )
  Sys.chmod(ds$backingfile,"0644")
  ds$attach(force=TRUE)
  expect_identical(ds[1, 1], 1, "OK read after chmod to readable and attach.")
  return(TRUE)
})

test_that("We can describe objects", {
    unlink(ds.data.file)
    mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
    ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
                                        # Getting
    expect_equal(nrow(ds),nrow(mat))
    expect_equal(ncol(ds),ncol(mat))
    expect_equal(dim(ds),dim(mat))
    expect_equal(length(ds),length(mat))
    expect_identical(dimnames(ds),dimnames(mat))
                                        # Setting
    new.dimnames = list(letters[4:6], LETTERS[4:6])
    dimnames(ds) = new.dimnames
    expect_identical(dimnames(ds), new.dimnames)

    unlink(ds.data.file)
    ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
    colnames(ds) = new.dimnames[[2]]
    expect_identical(colnames(ds), new.dimnames[[2]])

    unlink(ds.data.file)
    ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))
    rownames(ds) = new.dimnames[[1]]
    expect_identical(rownames(ds), new.dimnames[[1]])
})

test_that("We can reattach saved objects", {
  mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
  data.file = tempfile()
  old.ds = BigMatrix(mat,data.file,3,3,list(rownames,colnames))
  object.file = paste(data.file,".rds",sep="")
  saveRDS(old.ds,file=object.file)
  new.ds = readRDS(object.file)
  expect_equal( new.ds[,], mat )
  file.rename(old.ds$backingfile,file.path(tempdir(),"shoe2"))
  expect_error(old.ds$attach(force=TRUE))
})

test_that("We can get and set the data paths", {
  rownames = letters[1:3]
  colnames = LETTERS[1:3]
  mat = matrix(as.numeric(1:9),ncol=3,dimnames=list(rownames,colnames))
  int.mat = matrix(c(rep(1L,5),rep(2L,4)),ncol=3,dimnames=list(rownames,colnames))
  levels = c("AA","BB")
  char.mat = matrix()

  back.dir = tempdir()
  ds.data.file = file.path(back.dir,"bigmat","ds")
  unlink(ds.data.file)
  ds = BigMatrix(mat,ds.data.file,3,3,list(rownames,colnames))

  expect_equal(normalizePath(ds$backingfile),normalizePath(ds.data.file))
  new.data.file = tempfile()
  new.ds = BigMatrix(mat,new.data.file,3,3,list(rownames,colnames))
  even.newer.data.file = tempfile()
  expect_error({new.ds$backingfile = even.newer.data.file})
  file.copy(ds.data.file,even.newer.data.file)
  new.ds$backingfile = even.newer.data.file
  even.newer.data.file = normalizePath(even.newer.data.file)
  expect_equal(normalizePath(new.ds$backingfile),even.newer.data.file)
  expect_equal(new.ds[,], ds[,])
  expect_equal(basename(new.ds$backingfile), new.ds$.description$filename, "filename in description gets set by backingfile")
  desc = new.ds$description
  desc$filename = basename(new.data.file)
  new.ds$description = desc
  expect_equal(basename(new.ds$backingfile), new.ds$.description$filename, "filename in backingfile gets set by description")
})
