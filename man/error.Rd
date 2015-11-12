\name{error}
\alias{error}
\title{Forecast error measure}
\description{
Computes the forecast error measure.
}
\usage{
error(forecast, forecastbench, true, insampletrue, method = c("me", "mpe", "mae", 
 "mse", "sse", "rmse", "mdae", "mdse", "mape", "mdape", "smape", 
  "smdape", "rmspe", "rmdspe", "mrae", "mdrae", "gmrae", 
   "relmae", "relmse", "mase", "mdase", "rmsse"), giveall = FALSE)
}
\arguments{
  \item{forecast}{Out-of-sample forecasted values.}
  \item{forecastbench}{Forecasted values using a benchmark method, such as random walk.}
  \item{true}{Out-of-sample holdout values.}
  \item{insampletrue}{Insample values.}
  \item{method}{Method of forecast error measure.}
  \item{giveall}{If \code{giveall = TRUE}, all error measures are provided.}
}
\details{
\emph{\bold{Bias measure}}:

If \code{method = "me"}, the forecast error measure is mean error.

If \code{method = "mpe"}, the forecast error measure is mean percentage error.

\emph{\bold{Forecast accuracy error measure}}:

If \code{method = "mae"}, the forecast error measure is mean absolute error.

If \code{method = "mse"}, the forecast error measure is mean square error.

If \code{method = "sse"}, the forecast error measure is sum square error.

If \code{method = "rmse"}, the forecast error measure is root mean square error.

If \code{method = "mdae"}, the forecast error measure is median absolute error.

If \code{method = "mape"}, the forecast error measure is mean absolute percentage error.

If \code{method = "mdape"}, the forecast error measure is median absolute percentage error.

If \code{method = "rmspe"}, the forecast error measure is root mean square percentage error.

If \code{method = "rmdspe"}, the forecast error measure is root median square percentage error.

\emph{\bold{Forecast accuracy symmetric error measure}}:

If \code{method = "smape"}, the forecast error measure is symmetric mean absolute percentage error.

If \code{method = "smdape"}, the forecast error measure is symmetric median absolute percentage error.

\emph{\bold{Forecast accuracy relative error measure}}:

If \code{method = "mrae"}, the forecast error measure is mean relative absolute error.

If \code{method = "mdrae"}, the forecast error measure is median relative absolute error.

If \code{method = "gmrae"}, the forecast error measure is geometric mean relative absolute error.

If \code{method = "relmae"}, the forecast error measure is relative mean absolute error.

If \code{method = "relmse"}, the forecast error measure is relative mean square error.

\emph{\bold{Forecast accuracy scaled error measure}}:

If \code{method = "mase"}, the forecast error measure is mean absolute scaled error.

If \code{method = "mdase"}, the forecast error measure is median absolute scaled error.

If \code{method = "rmsse"}, the forecast error measure is root mean square scaled error.
}
\value{
A numeric value.
}
\references{
P. A. Thompson (1990) "An MSE statistic for comparing forecast accuracy across series", \emph{International Journal of Forecasting}, \bold{6}(2), 219-227.

C. Chatfield (1992) "A commentary on error measures", \emph{International Journal of Forecasting}, \bold{8}(1), 100-102.

S. Makridakis (1993) "Accuracy measures: theoretical and practical concerns", \emph{International Journal of Forecasting}, \bold{9}(4), 527-529.

R. J. Hyndman and A. Koehler (2006) "Another look at measures of forecast accuracy", \emph{International Journal of Forecasting}, \bold{22}(3), 443-473.
}
\author{Han Lin Shang}
\examples{
# Forecast error measures can be categorized into three groups: (1) scale-dependent, 
# (2) scale-independent but with possible zero denominator, 
# (3) scale-independent with non-zero denominator.
error(forecast = 1:2, true = 3:4, method = "mae")
error(forecast = 1:5, forecastbench = 6:10, true = 11:15, method = "mrae")
error(forecast = 1:5, forecastbench = 6:10, true = 11:15, insampletrue = 16:20, 
	giveall = TRUE)
}
\keyword{methods}


