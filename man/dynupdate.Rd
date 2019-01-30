\name{dynupdate}
\alias{dynupdate}
\title{Dynamic updates via BM, OLS, RR and PLS methods}
\description{ 
Four methods, namely block moving (BM), ordinary least squares (OLS) regression, 
ridge regression (RR), penalized least squares (PLS) regression, were proposed to 
address the problem of dynamic updating, when partial data in the most recent curve 
are observed. 
}
\usage{
dynupdate(data, newdata = NULL, holdoutdata, method = c("ts", "block", 
 "ols", "pls", "ridge"), fmethod = c("arima", "ar", "ets", "ets.na", 
  "rwdrift", "rw"), pcdmethod = c("classical", "M", "rapca"), 
   ngrid = max(1000, ncol(data$y)), order = 6, 
    robust_lambda = 2.33, lambda = 0.01, value = FALSE, 
     interval = FALSE, level = 80, 
      pimethod = c("parametric", "nonparametric"), B = 1000)
}
\arguments{
  \item{data}{An object of class \code{sfts}.}
  \item{newdata}{A data vector of newly arrived observations.}
  \item{holdoutdata}{A data vector of holdout sample to evaluate point forecast accuracy.}
  \item{method}{Forecasting methods. The latter four can dynamically update point forecasts.}
  \item{fmethod}{Univariate time series forecasting methods used in \code{method = "ts"} 
                 or \code{method = "block"}.}
  \item{pcdmethod}{Method to use for principal components decomposition. Possibilities are "M", "rapca" and "classical".}
  \item{ngrid}{Number of grid points to use in calculations. Set to maximum of 1000 and \code{ncol(data$y)}.}
  \item{order}{Number of principal components to fit.}
  \item{robust_lambda}{Tuning parameter in the two-step robust functional principal component analysis, when \code{pcdmethod = "M"}.}
  \item{lambda}{Penalty parameter used in \code{method = "pls"} or \code{method = "ridge"}.}
  \item{value}{When \code{value = TRUE}, returns forecasts or when \code{value = FALSE}, returns forecast errors.}
  \item{interval}{When \code{interval = TRUE}, produces distributional forecasts.}
  \item{level}{Nominal coverage probability.}
  \item{pimethod}{Parametric or nonparametric method to construct prediction intervals.}
  \item{B}{Number of bootstrap samples.}  
}
\details{
This function is designed to dynamically update point and interval forecasts, when partial data in the most recent curve are observed.

If \code{method = "classical"}, then standard functional principal component decomposition is used, as described by Ramsay and Dalzell (1991). 

If \code{method = "rapca"}, then the robust principal component algorithm of Hubert, Rousseeuw and Verboven (2002) is used. 

If \code{method = "M"}, then the hybrid algorithm of Hyndman and Ullah (2005) is used. 
}
\value{
  \item{forecasts}{An object of class \code{fts} containing the dynamic updated point forecasts.}
  \item{bootsamp}{An object of class \code{fts} containing the bootstrapped point forecasts, which are updated by the PLS method.}
  \item{low}{An object of class \code{fts} containing the lower bound of prediction intervals.}
  \item{up}{An object of class \code{fts} containing the upper bound of prediction intervals.}
}
\references{
J. O. Ramsay and C. J. Dalzell (1991) "Some tools for functional data analysis (with discussion)", \emph{Journal of the Royal Statistical Society: Series B}, \bold{53}(3), 539-572. 

M. Hubert and P. J. Rousseeuw and S. Verboven (2002) "A fast robust method for principal components with applications to chemometrics", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{60}(1-2), 101-111. 

R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

H. Shen and J. Z. Huang (2008) "Interday forecasting and intraday updating of call center arrivals", \emph{Manufacturing and Service Operations Management}, \bold{10}(3), 391-410.

H. Shen (2009) "On modeling and forecasting time series of curves", \emph{Technometrics}, \bold{51}(3), 227-238.

H. L. Shang and R. J. Hyndman (2011) "Nonparametric time series forecasting with dynamic updating", Mathematics and Computers in Simulation, \bold{81}(7), 1310-1324.

H. L. Shang (2013) "Functional time series approach for forecasting very short-term electricity demand", \emph{Journal of Applied Statistics}, \bold{40}(1), 152-168.

H. L. Shang (2017) "Forecasting Intraday S&P 500 Index Returns: A Functional Time Series Approach", \emph{Journal of Forecasting}, \bold{36}(7), 741-755.

H. L. Shang (2017) "Functional time series forecasting with dynamic updating: An application to intraday particulate matter concentration", \emph{Econometrics and Statistics}, \bold{1}, 184-200.
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}}
\examples{  
# ElNino is an object of sliced functional time series, constructed from a univariate time series. 
# When we observe some newly arrived information in the most recent time period, this function  
# allows us to update the point and interval forecasts for the remaining time period. 
dynupdate(data = ElNino_ERSST_region_1and2, newdata = ElNino_ERSST_region_1and2$y[1:4,69], 
	holdoutdata = ElNino_ERSST_region_1and2$y[5:12,57], method = "block", interval = FALSE) 
}
\keyword{models}

