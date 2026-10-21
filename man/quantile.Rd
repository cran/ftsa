\name{quantile}
\alias{quantile}
\title{Quantile}
\usage{
quantile(x, ...)
}
\arguments{
\item{x}{Numeric vector whose sample quantiles are wanted, or an object of a class for which a method has been defined.}
\item{...}{Arguments passed to specific methods.}
}
\description{Generic functions for quantile.}
\value{Refer to specific methods. For numeric vectors, see the \code{\link[stats]{quantile}} functions in the \code{stats} package.}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{quantile.fts}}}
\examples{
quantile(1:10)
quantile(x = ElNino_ERSST_region_1and2, probs = c(0.25, 0.5, 0.75))
}
\keyword{models}