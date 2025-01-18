\name{long_run_covariance_estimation}
\alias{long_run_covariance_estimation}
\title{
Estimating long-run covariance function for a functional time series
}
\description{
Bandwidth estimation in the long-run covariance function for a functional time series, using different types of kernel function
}
\usage{
long_run_covariance_estimation(dat, C0 = 3, H = 3)
}
\arguments{
  \item{dat}{A matrix of p by n, where p denotes the number of grid points and n denotes sample size}
  \item{C0}{A tuning parameter used in the adaptive bandwidth selection algorithm of Rice}
  \item{H}{A tuning parameter used in the adaptive bandwidth selection algorithm of Rice}
}
\value{
An estimated covariance function of size (p by p)
}
\references{
L. Horvath, G. Rice and S. Whipple (2016) Adaptive bandwidth selection in the long run covariance estimation of functional time series, \emph{Computational Statistics and Data Analysis}, \bold{100}, 676-693.

G. Rice and H. L. Shang (2017) A plug-in bandwidth selection procedure for long run covariance estimation with stationary functional time series, \emph{Journal of Time Series Analysis}, \bold{38}(4), 591-609.

D. Li, P. M. Robinson and H. L. Shang (2018) Long-range dependent curve time series, \emph{Journal of the American Statistical Association: Theory and Methods}, under revision.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[rainbow]{fts}}
}
\examples{
dum = long_run_covariance_estimation(dat = ElNino_OISST_region_1and2$y[,1:5])
}
\keyword{method}
