\name{BigMatrixFactor}
\alias{BigMatrixFactor}
\title{Create a new BigMatrixFactor}
\usage{
  BigMatrixFactor(x = NULL, backingfile, nrow, ncol,
    dimnames, levels)
}
\arguments{
  \item{x}{scalar numeric, NULL, matrix, or big.matrix.
  Optional data or big.matrix for new BigMatrix. A scalar
  numeric can be used to initalize the whole matrix.}

  \item{backingfile}{character, full path to the file that
  will contain the data matrix}

  \item{nrow}{integer, number of rows in the matrix we are
  about to create}

  \item{ncol}{integer, number of columns in the matrix we
  are about to create}

  \item{dimnames}{list, list(rownames,colnames), as for a
  typical matrix}

  \item{levels}{character, as for a typical factor}
}
\value{
  BigMatrixFactor
}
\description{
  Create a new BigMatrixFactor
}
\examples{
dnames = dimnames=list(letters[1:3],LETTERS[1:3])
x = matrix( sample( 1:3, 9, replace=TRUE), ncol=3, dimnames=dnames)
ds = BigMatrixFactor(x,tempfile(),levels=c("AA","AB","BB"))
ds = BigMatrixFactor(backingfile=tempfile(),nrow=3,ncol=3,dimnames=dnames,levels=c("AA","AB","BB"))
}

