\name{H_estimate_LW_NFTS}
\alias{H_estimate_LW_NFTS}
\title{
Estimation of the long-memory parameter in nonstationary functional time series
}
\description{
Estimation of the long-memory parameter via local Whittle estimator in nonstationary functional time series
}
\usage{
H_estimate_LW_NFTS(fun_dat, band_const, tuning_const, delta_val, m_const, 
est_method = c("diff_add", "ELW", "ELW_2S"), choice_score = c("first", "average"), 
taper_choice = c("none", "Velasco", "HC"), 
ELW_mean_est = c("mean", "init", "weighted", "none"))
}
\arguments{
  \item{fun_dat}{A functional data set}
  \item{band_const}{Tuning parameter, customarily set to 1}
  \item{tuning_const}{Tuning parameter, customarily set to 0.4}
  \item{delta_val}{Tuning parameter, customarily set to 0.01}
  \item{m_const}{Tuning parameter, customarily set to 0.65}
  \item{est_method}{Estimation method from a set of principal component scores}
  \item{choice_score}{The leading set of principal component scores or a norm}
  \item{taper_choice}{Tapering parameter, customarily none}
  \item{ELW_mean_est}{When the estimation method is ELW, it is the mean term of the ELW}
}
\value{
\item{C}{Estimated long-run covariance}
\item{score}{Estimated principal component scores}
\item{H_est}{Hurst parameter estimate}
\item{d_est}{Long-memory parameter estimate}
\item{alpha_est}{Long-memory parameter estimate}
\item{ncomp}{Number of retained components}
\item{ER}{Eigenvalue ratio}
\item{diff_order}{Order of difference}
\item{eigen_function}{Principal components}
}
\references{
D. Li, P. M. Robinson, H. L. Shang (2020) Long-range dependent curve time series, \emph{Journal of the American Statistical Association}, \bold{115}(530), 957-971. 

D. Li, P. M. Robinson, H. L. Shang (2021) Local Whittle estimation of long-range dependence for functional time series, \emph{Journal of Time Series Analysis}, \bold{42}(5-6), 685-695.

D. Li, P. M. Robinson, H. L. Shang (2023) Nonstationary fractionally integrated functional time series, \emph{Bernoulli}, \bold{29}(2), 1505-1526.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{H_estimate_RS}}
}
\examples{
sim_data_set = sim_FARMA(n = 50, seed_number = 123, process = "FAR")
LW_d_est = H_estimate_LW_NFTS(fun_dat=sim_data_set, band_const = 1, tuning_const = 0.4, 
	delta_val = 0.01, m_const = 0.65, est_method = "diff_add", choice_score = "first", 
	taper_choice = "none")
}
\keyword{methods}
