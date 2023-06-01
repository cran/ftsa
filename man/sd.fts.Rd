\name{sd.fts}
\alias{sd.fts}
\title{Standard deviation functions for functional time series}
\description{
Computes standard deviation of functional time series at each variable.
}
\usage{
\method{sd}{fts}(x, method = c("coordinate", "FM", "mode", "RP", "RPD", "radius"), 
 trim = 0.25, alpha, weight,...)
}
\arguments{
  \item{x}{An object of class \code{fts}.}
  \item{method}{Method for computing median.}
  \item{trim}{Percentage of trimming.}
  \item{alpha}{Tuning parameter when \code{method="radius"}.}
  \item{weight}{Hard thresholding or soft thresholding.}
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

If \code{method = "RPD"}, it computes the standard deviation functions of trimmed functional data ordered by random projection with
derivative depth.

If \code{method = "radius"}, it computes the standard deviation function of trimmed functional data ordered by the notion of alpha-radius.
}
\value{
A list containing \code{x} = variables and \code{y} = standard deviation rates.
}
\references{
O. Hossjer and C. Croux (1995) "Generalized univariate signed rank statistics for testing and estimating a multivariate location parameter", \emph{Nonparametric Statistics}, \bold{4}(3), 293-308.

A. Cuevas and M. Febrero and R. Fraiman (2006) "On the use of bootstrap for estimating functions with functional data", \emph{Computational Statistics and Data Analysis}, \bold{51}(2), 1063-1074.

A. Cuevas and M. Febrero and R. Fraiman (2007), "Robust estimation and classification for functional data via projection-based depth notions", \emph{Computational Statistics}, \bold{22}(3), 481-496.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2007) "A functional analysis of NOx levels: location and scale estimation and outlier detection", \emph{Computational Statistics}, \bold{22}(3), 411-427.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2010) "Measures of influence for the functional linear model with scalar response", \emph{Journal of Multivariate Analysis}, \bold{101}(2), 327-339. 

J. A. Cuesta-Albertos and A. Nieto-Reyes (2010) "Functional classification and the random Tukey depth. Practical issues", Combining Soft Computing and Statistical Methods in Data Analysis, \emph{Advances in Intelligent and Soft Computing}, \bold{77}, 123-130.

D. Gervini (2012) "Outlier detection and trimmed estimation in general functional spaces", \emph{Statistica Sinica}, \bold{22}(4), 1639-1660.
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{median.fts}}, \code{\link[ftsa]{var.fts}}, 
\code{\link[ftsa]{quantile.fts}}}
\examples{
# Fraiman-Muniz depth was arguably the oldest functional depth.	
sd(x = ElNino_ERSST_region_1and2, method = "FM")
sd(x = ElNino_ERSST_region_1and2, method = "coordinate")
sd(x = ElNino_ERSST_region_1and2, method = "mode")
sd(x = ElNino_ERSST_region_1and2, method = "RP")
sd(x = ElNino_ERSST_region_1and2, method = "RPD")
sd(x = ElNino_ERSST_region_1and2, method = "radius", 
	alpha = 0.5, weight = "hard")
sd(x = ElNino_ERSST_region_1and2, method = "radius", 
	alpha = 0.5, weight = "soft")
}
\keyword{methods}

