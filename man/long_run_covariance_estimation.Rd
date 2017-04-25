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
L. Horvath, G. Rice and S. Whipple (2016) Adaptive bandwidth selection in the long run covariance estimation of functional time series, Computational Statistics and Data Analysis, 100, 676-693.

G. Rice and H. L. Shang (2016) A plug-in bandwidth selection procedure for long run covariance estimation with stationary functional time series, Journal of Time Series Analysis, in press.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{fts}}
}
\examples{
dum = long_run_covariance_estimation(dat = ElNino$y)
}
\keyword{method}
