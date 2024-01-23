\name{facf}
\alias{facf}
\title{
Functional autocorrelation function
}
\description{
Compute functional autocorrelation function at various lags
}
\usage{
facf(fun_data, lag_value_range = seq(0, 20, by = 1))
}
\arguments{
  \item{fun_data}{A data matrix of dimension (n by p), where n denotes sample size; and p denotes dimensionality}
  \item{lag_value_range}{Lag value}
}
\details{
The autocovariance at lag \eqn{i} is estimated by the function \eqn{\widehat{\gamma}_i(t,s)}, a functional analog of the autocorrelation is defined as
\deqn{\widehat{\rho}_i = \frac{\|\widehat{\gamma}_i\|}{\int \widehat{\gamma}_0(t,t)dt}.}{\widehat{\rho}_i = \frac{\|\widehat{\gamma}_i\|}{\int \widehat{\gamma}_0(t,t)dt}.}
}
\value{
A vector of functional autocorrelation function at various lags
}
\references{
L. Horv\'ath, G. Rice and S. Whipple (2016) Adaptive bandwidth selection in the long run covariance estimator of functional time series, \emph{Computational Statistics and Data Analysis}, \bold{100}, 676-693.
}
\author{
Han Lin Shang
}
\examples{
facf_value = facf(fun_data = t(ElNino_ERSST_region_1and2$y))
}
\keyword{methods}
