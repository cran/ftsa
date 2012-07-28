\name{sd.fts}
\alias{sd.fts}
\title{Standard deviation functions for functional time series}
\description{
Computes standard deviation of functional time series at each variable.
}
\usage{
\method{sd}{fts}(x, method = c("coordinate", "FM", "mode", "RP", "RPD"), 
 trim = 0.25,...)
}
\arguments{
  \item{x}{An object of class \code{fts}.}
  \item{method}{Method for computing median.}
  \item{trim}{Percentage of trimming.}
  \item{...}{Other arguments.}
}
\details{
If \code{method = "coordinate"}, it computes coordinate-wise standard deviation functions.

If \code{method = "FM"}, it computes the standard deviation functions of trimmed functional data ordered by the functional depth of 
Fraiman and Muniz (2001).

If \code{method = "mode"}, it computes the standard deviation functions of trimmed functional data ordered by \eqn{h}-modal functional 
depth.

If \code{method = "RP"}, it computes the standard deviation functions of trimmed functional data ordered by random projection 
depth.

If \code{method = "RPD"}, it computes the standard deviation functions of trimmed functional data ordered by random projection 
derivative depth.
}
\value{
A list containing \code{x} = variables and \code{y} = standard deviation rates.
}
\references{
O. Hossjer and C. Croux (1995) "Generalized univariate signed rank statistics for testing and estimating a multivariate location parameter", \emph{Nonparametric Statistics}, \bold{4}(3), 293-308.

A. Cuevas and M. Febrero and R. Fraiman (2006) "On the use of bootstrap for estimating functions with functional data", \emph{Computational Statistics \& Data Analysis}, \bold{51}(2), 1063-1074.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{median.fts}}, \code{\link[ftsa]{var.fts}}, 
\code{\link[ftsa]{quantile.fts}}}
\examples{
sd(x = ElNino, method = "FM")
}
\keyword{methods}

