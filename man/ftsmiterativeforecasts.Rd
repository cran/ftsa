\name{ftsmiterativeforecasts}
\alias{ftsmiterativeforecasts}
\title{Forecast functional time series}
\description{The coefficients from the fitted object are forecasted
using either an ARIMA model (\code{method = "arima"}), an AR model (\code{method = "ar"}), 
an exponential smoothing method (\code{method = "ets"}), a linear exponential smoothing 
method allowing missing values (\code{method = "ets.na"}), or a random walk with drift model 
(\code{method = "rwdrift"}). The forecast coefficients are then multiplied by the principal 
components to obtain a forecast curve.
}
\usage{
ftsmiterativeforecasts(object, components, iteration = 20)
}
\arguments{
  \item{object}{An object of class \code{fts}.}
  \item{components}{Number of principal components.}
  \item{iteration}{Number of iterative one-step-ahead forecasts.}
}
\details{
1. Obtain a smooth curve \eqn{f_t(x)} for each \eqn{t} using a nonparametric smoothing technique.

2. Decompose the smooth curves via a functional principal component analysis.

3. Fit a univariate time series model to each of the principal component scores.

4. Forecast the principal component scores using the fitted time series models.

5. Multiply the forecast principal component scores by fixed principal components to obtain forecasts of \eqn{f_{n+h}(x)}.

6. The estimated variances of the error terms (smoothing error and model residual error) are used to compute prediction intervals for the forecasts.
}
\value{
List with the following components:
\item{mean}{An object of class \code{fts} containing point forecasts.}
\item{lower}{An object of class \code{fts} containing lower bound for prediction intervals.}
\item{upper}{An object of class \code{fts} containing upper bound for prediction intervals.}
\item{fitted}{An object of class \code{fts} of one-step-ahead forecasts for historical data.}
\item{error}{An object of class \code{fts} of one-step-ahead errors for historical data.}
\item{coeff}{List of objects of type \code{forecast} containing the coefficients and their forecasts.}
\item{coeff.error}{One-step-ahead forecast errors for each of the coefficients.}
\item{var}{List containing the various components of variance: model, error, mean, total and coeff.}
\item{model}{Fitted \code{\link[ftsa]{ftsm}} model.}
\item{bootsamp}{An array of \eqn{dim = c(p, B, h)} containing the bootstrapped point forecasts. 
                \eqn{p} is the number of variables. \eqn{B} is the number of bootstrap samples.
                \eqn{h} is the forecast horizon.}
}
\references{
H. Booth and R. J. Hyndman and L. Tickle and P. D. Jong (2006) "Lee-Carter mortality forecasting: A multi-country comparison of variants and extensions", \emph{Demographic Research}, \bold{15}, 289-310.

B. Erbas and R. J. Hyndman and D. M. Gertig (2007) "Forecasting age-specific breast cancer mortality using functional data model", \emph{Statistics in Medicine}, \bold{26}(2), 458-470.

R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series" (with discussion), \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{
Han Lin Shang
}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}}
\examples{
# Iterative one-step-ahead forecasts via functional principal component analysis.	
ftsmiterativeforecasts(object = Australiasmoothfertility, components = 2, iteration = 5)
}
\keyword{models}
