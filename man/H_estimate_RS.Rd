\name{H_estimate_RS}
\alias{H_estimate_RS}
\title{
Estimation of long-memory parameter
}
\description{
R/S estimate of the long-memory functional time series
}
\usage{
H_estimate_RS(fun_dat, band_const = 1, choice_score = c("first", "average"),
    est_method = c("classical", "modified"))
}
\arguments{
  \item{fun_dat}{Functional data object on a set of discretized grid points}
  \item{band_const}{Constant multiplier}
  \item{choice_score}{First set of principal component scores or the norm of multiple principal component scores}
  \item{est_method}{Classical or modified R/S estimator}
}
\value{
\item{C}{Estimated covariance}
\item{score}{Estimated principal component score}
\item{H_est}{Estimated Hurst exponent}
\item{d_est}{Estimated memory parameter}
\item{alpha_est}{Estimated alpha parameter}
\item{ncomp}{Number of retained components}
\item{prop_eigen_value}{Proportions of eigenvalues}
\item{eigen_function}{Retained principal components}
}
\references{
D. Li, P. M. Robinson, H. L. Shang (2020) Long-range dependent curve time series. Journal of the American Statistical Association, 115(530), 957-971.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{sim_FARMA}}
}
\examples{
sim_data_set = sim_FARMA(n = 50, seed_number = 123, process = "FAR")
RS_est = H_estimate_RS(sim_data_set, choice_score = "first", est_method = "classical")
}
\keyword{time series}
