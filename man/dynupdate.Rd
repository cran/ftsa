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
  "rwdrift", "rw"), error = c("mse", "mae", "mape"), order = 6,
   lambda = 0.01, value = FALSE, interval = FALSE, 
    level = 80, B = 1000)
}
\arguments{
  \item{data}{An object of class \code{sfts}.}
  \item{newdata}{A data vector of newly arrived observations.}
  \item{holdoutdata}{A data vector of holdout sample to evaluate point forecast accuracy.}
  \item{method}{Forecasting methods. The latter four can dynamically update point forecasts.}
  \item{fmethod}{Univariate time series forecasting methods used in \code{method = "ts"} 
                 or \code{method = "block"}.}
  \item{error}{Forecast error measures.}
  \item{order}{Number of principal components to fit.}
  \item{lambda}{Penalty parameter used in \code{method = "pls"}.}
  \item{value}{When \code{value = TRUE}, returns forecasts or when \code{value = FALSE}, returns forecast errors.}
  \item{interval}{When \code{interval = TRUE}, produces distributional forecasts.}
  \item{level}{Nominal coverage probability.}
  \item{B}{Number of bootstrap samples.}  
}
\details{
This function is designed to dynamically update point and distributional forecasts, when partial data in the most recent curve are observed.
}
\value{
  \item{forecasts}{An object of class \code{fts} containing the dynamic updated point forecasts.}
  \item{bootsamp}{An object of class \code{fts} containing the bootstrapped point forecasts, which are updated by the PLS method.}
  \item{low}{An object of class \code{fts} containing the lower bound of prediction intervals.}
  \item{up}{An object of class \code{fts} containing the upper bound of prediction intervals.}
}
\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics \& Data Analysis}, \bold{51}(10), 4942-4956.

H. Shen and J. Z. Huang (2008) "Interday forecasting and intraday updating of call center arrivals", \emph{Manufacturing \& Service Operations Management}, \bold{10}(3), 391-410.

H. Shen (2009) "On modeling and forecasting time series of curves", \emph{Technometrics}, \bold{51}(3), 227-238.

H. Shang and R. J. Hyndman (2009) "Nonparametric time series forecasting with dynamic updating", In R. S. Anderssen, R. D. Braddock and L.T.H. Newham (eds), 18th World IMACS Congress and MODSIM09 International Congress on Modelling and Simulation. Modelling and Simulation Society of Australia and New Zealand and International Association for Mathematics and Computers in Simulation, July 2009, pp. 1552-1558. ISBN: 978-0-9758400-7-8. 
\url{http://www.mssanz.org.au/modsim09/D11/shang.pdf}
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}}
\examples{
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "ts")
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "block")
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "ols")
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "pls")
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "ridge")
dynupdate(data = ElNino, newdata = ElNino$y[1:4,57], 
          holdoutdata = ElNino$y[5:12,57], method = "block", interval = TRUE)
}
\keyword{models}

