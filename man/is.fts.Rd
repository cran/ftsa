\name{is.fts}
\alias{is.fts}
\title{Test for functional time series}
\description{
Tests whether an object is of class \code{fts}. 
}
\usage{
is.fts(x)
}
\arguments{
  \item{x}{Arbitrary R object.}
}
\author{Rob J Hyndman}
\examples{
# check if ElNino is the class of the functional time series.	
is.fts(x = ElNino)
}
\keyword{ts}


