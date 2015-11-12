\name{farforecast}
\alias{farforecast}
\title{
Functional data forecasting through functional principal component regression
}
\description{
The coefficients from the fitted object are forecasted using a multivariate time-series forecasting method. 
The forecast coefficients are then multiplied by the functional principal components to obtain a forecast curve.
}
\usage{
farforecast(object, h = 10, var_type = "const", level = 80, PI = FALSE)
}
\arguments{
  \item{object}{An object of \code{\link{ftsm}}.}
  \item{h}{Forecast horizon.}
  \item{var_type}{Type of multivariate time series forecasting method; see \code{\link[vars]{VAR}} for details.}
  \item{level}{Nominal coverage probability.}
  \item{PI}{When \code{PI = TRUE}, a prediction interval will be given along with the point forecast.}
}
\details{
1. Decompose the smooth curves via a functional principal component analysis (FPCA).

2. Fit a multivariate time-series model to the principal component score matrix.

3. Forecast the principal component scores using the fitted multivariate time-series models. The order of VAR is selected optimally via an information criterion.

4. Multiply the forecast principal component scores by estimated principal components to obtain forecasts of \eqn{f_{n+h}(x)}.

5. Prediction intervals are constructed by taking quantiles of the one-step-ahead forecast errors. 
}
\value{
\item{point_fore}{Point forecast}
\item{PI_lb}{Lower bound of a prediction interval}
\item{PI_ub}{Upper bound of a prediction interval}
}
\references{
A. Aue, D. D. Norinho and S. Hormann (2015) "On the prediction of stationary functional time series", \emph{Journal of the American Statistical Association}, \bold{110}(509), 378-392.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{forecast.ftsm}}, \code{\link{forecastfplsr}}
}
\examples{
ElNino_subset = extract(data = ElNino, direction = "time",timeorder = 1967:2006)
ex = farforecast(ftsm(y = ElNino_subset), PI = TRUE)
}
\keyword{models}