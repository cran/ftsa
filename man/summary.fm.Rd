\name{summary.fm}
\alias{summary.fm}
\title{Summary for functional time series model}
\description{
Summarizes a basis function model fitted to a functional time series. It returns various measures of goodness-of-fit.
}
\usage{
\method{summary}{fm}(object, ...)
}
\arguments{
\item{object}{Output from \code{\link[ftsa]{ftsm}} or \code{\link[ftsa]{fplsr}}.}
\item{...}{Other arguments.}
}
\value{None.}

\author{Rob J Hyndman}

\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.fmres}}}
\examples{
summary(object = ftsm(ElNino))
}
\keyword{models}
