% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collections.R
\name{updateBackingfiles}
\alias{updateBackingfiles}
\title{Update directory for BigMatrix objects in a collection to new location}
\usage{
updateBackingfiles(x, dir)
}
\arguments{
\item{x}{list, SimpleList, environment or something with names and [[ methods}

\item{dir}{character, path to directory holding all BigMatrix files}
}
\value{
x param, with modified BigMatrix objects.
}
\description{
Update directory for BigMatrix objects in a collection to new location. Assumes files have already
been moved on the filesystem. Assumes names of description and data files are the same. The
collection can contain a mix of BigMatrix objects and other types. The other types will be
not be touched and will be returned as they are.
}
\details{
If you have renamed specific backing files, you will want to update the
backingfile field of the relevant BigMatrix objects.
}
\examples{
\dontrun{ list = updateBackingfiles(list, "/new/path/to/bigmat/dir") }
\dontrun{ assays(se) = updateBackingfiles(assays(se), "/new/path/to/bigmat/dir") }
\dontrun{ assayData(eset) = updateBackingfiles(assayData(eset), "/new/path/to/bigmat/dir") }
}
