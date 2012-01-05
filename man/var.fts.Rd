\name{var.fts}
\alias{var.fts}
\title{Variance functions for functional time series}
\description{
Computes variance functions of functional time series at each variable.
}
\usage{
\method{var}{fts}(x, method = c("coordinate", "FM", "mode", "RP", "RPD"), 
 trim = 0.25,...)
}
\arguments{
  \item{x}{An object of class \code{fts}.}
  \item{method}{Method for computing median.}
  \item{trim}{Percentage of trimming.}
  \item{...}{Other arguments.}
}
\details{
If \code{method = "coordinate"}, it computes coordinate-wise variance.

If \code{method = "FM"}, it computes the variance of trimmed functional data ordered by the functional depth of Fraiman and Muniz (2001).

If \code{method = "mode"}, it computes the variance of trimmed functional data ordered by \eqn{h}-modal functional depth.

If \code{method = "RP"}, it computes the variance of trimmed functional data ordered by random projection depth.

If \code{method = "RPD"}, it computes the variance of trimmed functional data ordered by random projection derivative depth.
}
\value{
A list containing \code{x} = variables and \code{y} = variance rates.
}
\references{
O. Hossjer and C. Croux (1995) "Generalized univariate signed rank statistics for testing and estimating a multivariate location parameter", \emph{Nonparametric Statistics}, \bold{4}(3), 293-308.

A. Cuevas and M. Febrero and R. Fraiman (2006) "On the use of bootstrap for estimating functions with functional data", \emph{Computational Statistics \& Data Analysis}, \bold{51}(2), 1063-1074.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{median.fts}}, \code{\link[ftsa]{sd.fts}}, \code{\link[ftsa]{quantile.fts}}}
\examples{
var(x = ElNino, method = "coordinate")
var(x = ElNino, method = "FM")
var(x = ElNino, method = "mode")
var(x = ElNino, method = "RP")
var(x = ElNino, method = "RPD")
}
\keyword{methods}

