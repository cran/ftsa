\name{forecast.hdfpca}
\alias{forecast.hdfpca}
\title{
Forecasting via a high-dimensional functional principal component regression
}
\description{
Forecast high-dimensional functional principal component model.
}
\usage{
\method{forecast}{hdfpca}(object, h = 3, level = 80, B = 50, ...)
}
\arguments{
  \item{object}{An object of class 'hdfpca'}
  \item{h}{Forecast horizon}
  \item{level}{Prediction interval level, the default is 80 percent}
  \item{B}{Number of bootstrap replications}
  \item{...}{Other arguments passed to forecast routine.}
}
\details{
The low-dimensional factors are forecasted with autoregressive integrated moving average (ARIMA) models separately. The forecast functions are then calculated using the forecast factors. Bootstrap prediction intervals are constructed by resampling from the forecast residuals of the ARIMA models.
}
\value{
\item{forecast}{A list containing the h-step-ahead forecast functions for each population}
\item{upper}{Upper confidence bound for each population}
\item{lower}{Lower confidence bound for each population}
}
\references{
Y. Gao, H. L. Shang and Y. Yang (2018) High-dimensional functional time series forecasting: An application to age-specific mortality rates, \emph{Journal of Multivariate Analysis}, \bold{forthcoming}.
}
\author{
Y. Gao and H. L. Shang
}
\seealso{
\code{\link[ftsa]{hdfpca}}, \code{\link[ftsa]{hd_data}} 
}
\examples{
\dontrun{
hd_model = hdfpca(hd_data, order = 2, r = 2)
hd_model_fore = forecast.hdfpca(object = hd_model, h = 1)
}
}
\keyword{models}