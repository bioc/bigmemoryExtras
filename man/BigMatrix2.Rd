\name{BigMatrix2}
\alias{BigMatrix2}
\title{Create a new BigMatrix2}
\usage{
  BigMatrix2(x = NULL, backingfile, nrow, ncol,
    dimnames = NULL, type = "double")
}
\arguments{
  \item{x}{scalar numeric, NULL, matrix, or big.matrix.
  Optional data or big.matrix for new BigMatrix2. A scalar
  numeric can be used to initalize the whole matrix.}

  \item{backingfile}{character, full path to the file that
  will contain the data matrix}

  \item{nrow}{integer, number of rows in the matrix we are
  about to create}

  \item{ncol}{integer, number of columns in the matrix we
  are about to create}

  \item{dimnames}{list, list(rownames,colnames), as for a
  typical matrix}

  \item{type}{character type of big.matrix (double,
  integer, char)}
}
\value{
  BigMatrix2
}
\description{
  Create a new BigMatrix2
}
\examples{
dnames = dimnames=list(letters[1:3],LETTERS[1:3])
x = matrix(1:9,ncol=3,dimnames=dnames)
ds = BigMatrix2(x,tempfile())
ds = BigMatrix2(backingfile=tempfile(),nrow=3,ncol=3,dimnames=dnames)
}

