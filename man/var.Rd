\name{var}
\alias{var}
\alias{var.default}
\title{Variance}
\usage{
var(...)
}
\arguments{
\item{...}{Arguments passed to specific methods.}
}
\description{Generic functions for variance.}
\details{The \code{\link[stats]{var}} functions in the \code{stats} package are replaced by \code{var.default}.}
\value{Refer to specific methods. For numeric vectors, see the \code{\link[stats]{var}} functions in the \code{stats} 
package.}
\author{Rob J Hyndman and Han Lin Shang}
\seealso{\code{\link[ftsa]{var.fts}}}
\examples{
var(1:10)
var(x = ElNino_ERSST_region_1and2, method = "coordinate")
}
\keyword{models}
