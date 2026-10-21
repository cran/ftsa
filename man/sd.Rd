\name{sd}
\alias{sd}
\alias{sd.default}
\title{Standard deviation}
\usage{
sd(...)
}

\arguments{
\item{...}{Arguments passed to specific methods.}
}

\description{Generic functions for standard deviation.}

\details{The \code{\link[stats]{sd}} functions in the \code{stats} package are replaced by \code{sd.default}.}

\value{Refer to specific methods. For numeric vectors, see the \code{\link[stats]{sd}} functions in the \code{stats}
package.}

\author{Han Lin Shang}

\seealso{\code{\link[ftsa]{sd.fts}}}
\examples{
# the generic dispatches to sd.default for a numeric vector
sd(1:10)
# and to sd.fts for a functional time series
sd(x = ElNino_ERSST_region_1and2, method = "coordinate")
}
\keyword{models}
