\name{dynamic_FLR}
\alias{dynamic_FLR}
\title{
Dynamic updates via functional linear regression
}
\description{
A functional linear regression is used to address the problem of dynamic updating, when partial data in the most recent curve are observed.
}
\usage{
dynamic_FLR(dat, newdata, holdoutdata, order_k_percent = 0.9, order_m_percent = 0.9, 
    pcd_method = c("classical", "M"), robust_lambda = 2.33, bootrep = 100, 
    	pointfore, level = 80)
}
\arguments{
  \item{dat}{An object of class \code{sfts}.}
  \item{newdata}{A data vector of newly arrived observations.}
  \item{holdoutdata}{A data vector of holdout sample to evaluate point forecast accuracy.}
  \item{order_k_percent}{Select the number of components that explains at least 90 percent of the total variation.}
  \item{order_m_percent}{Select the number of components that explains at least 90 percent of the total variation.}
  \item{pcd_method}{Method to use for principal components decomposition. Possibilities are "M", "rapca" and "classical".}
  \item{robust_lambda}{Tuning parameter in the two-step robust functional principal component analysis, when \code{pcdmethod = "M"}.}
  \item{bootrep}{Number of bootstrap samples.}
  \item{pointfore}{If \code{pointfore = TRUE}, point forecasts are produced.}
  \item{level}{Nominal coverage probability.}
}
\details{
This function is designed to dynamically update point and interval forecasts, when partial data in the most recent curve are observed.
}
\value{
\item{update_forecast}{Updated forecasts.}
\item{holdoutdata}{Holdout sample.}
\item{err}{Forecast errors.}
\item{order_k}{Number of principal components in the first block of functions.}
\item{order_m}{Number of principal components in the second block of functions.}
\item{update_comb}{Bootstrapped forecasts for the dynamically updating time period.}
\item{update_comb_lb_ub}{By taking corresponding quantiles, obtain lower and upper prediction bounds.}
\item{err_boot}{Bootstrapped in-sample forecast error for the dynamically updating time period.}
}
\references{
H. Shen and J. Z. Huang (2008) "Interday forecasting and intraday updating of call center arrivals", \emph{Manufacturing and Service Operations Management}, \bold{10}(3), 391-410.

H. Shen (2009) "On modeling and forecasting time series of curves", \emph{Technometrics}, \bold{51}(3), 227-238.

H. L. Shang and R. J. Hyndman (2011) "Nonparametric time series forecasting with dynamic updating", Mathematics and Computers in Simulation, \bold{81}(7), 1310-1324.

J-M. Chiou (2012) "Dynamical functional prediction and classification with application to traffic flow prediction", \emph{Annals of Applied Statistics}, \bold{6}(4), 1588-1614.

H. L. Shang (2013) "Functional time series approach for forecasting very short-term electricity demand", \emph{Journal of Applied Statistics}, \bold{40}(1), 152-168.

H. L. Shang (2015) "Forecasting Intraday S&P 500 Index Returns: A Functional Time Series Approach", \emph{Journal of Forecasting}, \bold{36}(7), 741-755.

H. L. Shang (2017) "Functional time series forecasting with dynamic updating: An application to intraday particulate matter concentration", \emph{Econometrics and Statistics}, \bold{1}, 184-200.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{dynupdate}}
}
\examples{
dynamic_FLR_point = dynamic_FLR(dat = ElNino_ERSST_region_1and2$y[,1:68], 
	newdata = ElNino_ERSST_region_1and2$y[1:4,69], 
	holdoutdata = ElNino_ERSST_region_1and2$y[5:12,69], pointfore = TRUE)

dynamic_FLR_interval = dynamic_FLR(dat = ElNino_ERSST_region_1and2$y[,1:68], 
	newdata = ElNino_ERSST_region_1and2$y[1:4,69], 
	holdoutdata = ElNino_ERSST_region_1and2$y[5:12,69], pointfore = FALSE)
}
\keyword{models}
