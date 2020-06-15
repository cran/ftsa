\name{diff.fts}
\alias{diff.fts}
\title{Differences of a functional time series}
\description{
Computes differences of a \code{fts} object at each variable.
}
\usage{
\method{diff}{fts}(x, lag = 1, differences = 1, ...)
}
\arguments{
\item{x}{An object of class \code{fts}.}
\item{lag}{An integer indicating which lag to use.}
\item{differences}{An integer indicating the order of the difference.}
\item{...}{Other arguments.}
}
\value{An object of class \code{fts}.}

\author{Rob J Hyndman}
\examples{
# ElNino is an object of sliced functional time series.
# Differencing is sometimes used to achieve stationarity.	
diff(x = ElNino_ERSST_region_1and2)
}
\keyword{ts}
