\name{quantile.fts}
\alias{quantile.fts}
\title{Quantile functions for functional time series}
\description{
Computes quantiles of functional time series at each variable.
}
\usage{
\method{quantile}{fts}(x, probs, ...)
}
\arguments{
\item{x}{An object of class \code{fts}.}
\item{probs}{Quantile percentages.}
\item{...}{Other arguments.}
}
\value{
Return quantiles for each variable.
}
\author{Han Lin Shang}
\seealso{
\code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{median.fts}}, \code{\link[ftsa]{var.fts}}, \code{\link[ftsa]{sd.fts}}
}
\examples{
quantile(x = ElNino_ERSST_region_1and2)
}
\keyword{methods}

