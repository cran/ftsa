\name{extract}
\alias{extract}
\title{Extract variables or observations}
\description{
Creates subsets of a \code{fts} object.
}
\usage{
extract(data, direction = c("time", "x"), timeorder, xorder)
}
\arguments{
   \item{data}{An object of \code{fts}.}
   \item{direction}{In time direction or x variable direction?}
   \item{timeorder}{Indexes of time order.}
   \item{xorder}{Indexes of x variable order.}
}
\value{
When \code{xorder} is specified, it returns a \code{fts} object with same argument as 
data but with a subset of \code{x} variables.\cr

When \code{timeorder} is specified, it returns a \code{fts} object with same argument 
as data but with a subset of \code{time} variables.
}
\author{Han Lin Shang}
\examples{
# ElNino is an object of class sliced functional time series.
# This function truncates the data series rowwise or columnwise.	
extract(data = ElNino_ERSST_region_1and2, direction = "time", 
	timeorder = 1980:2006) # Last 27 curves
extract(data = ElNino_ERSST_region_1and2, direction = "x", 
	xorder = 1:8) # First 8 x variables
}
\keyword{models}

