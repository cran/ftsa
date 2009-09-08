\name{median.fts}
\alias{median.fts}
\title{Median functions for functional time series}
\description{
Computes median of functional time series at each variable.
}
\usage{
median.fts(x, method = c("hossjercroux", "coordinate", "FM", "mode", 
 "RP", "RPD"), ...)
}
\arguments{
  \item{x}{An object of class \code{fts}.}
  \item{method}{Method for computing median.}
  \item{\dots}{Other arguments.}
}
\details{
If \code{method = "coordinate"}, it computes a coordinate-wise median.

If \code{method = "hossjercroux"}, it computes the L1-median using the Hossjer-Croux algorithm.

If \code{method = "FM"}, it computes the median of trimmed functional data ordered by the functional depth of Fraiman and Muniz (2001).

If \code{method = "mode"}, it computes the median of trimmed functional data ordered by \eqn{h}-modal functional depth.

If \code{method = "RP"}, it computes the median of trimmed functional data ordered by random projection depth.

If \code{method = "RPD"}, it computes the median of trimmed functional data ordered by random projection derivative depth.
}
\value{
A list containing \code{x} = variables and \code{y} = median rates.
}
\references{
O. Hossjer and C. Croux (1995) "Generalized univariate signed rank statistics for testing and estimating a multivariate location parameter", \emph{Nonparametric Statistics}, \bold{4}(3), 293-308.

A. Cuevas and M. Febrero and R. Fraiman (2006) "On the use of bootstrap for estimating functions with functional data", \emph{Computational Statistics \& Data Analysis}, \bold{51}(2), 1063-1074.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.
}
\author{Rob J Hyndman, Han Lin Shang}
\seealso{\code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{var.fts}}, \code{\link[ftsa]{sd.fts}}, \code{\link[ftsa]{quantile.fts}}}
\examples{
median(x = ausmale, method = "hossjercroux")
median(x = ausmale, method = "coordinate")
median(x = ausmale, method = "FM")
median(x = ausmale, method = "mode")
median(x = ausmale, method = "RP")
median(x = ausmale, method = "RPD")
}
\keyword{methods}

