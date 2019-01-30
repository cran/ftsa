\name{forecast.ftsm}
\alias{forecast.ftsm}
\title{Forecast functional time series}
\description{The coefficients from the fitted object are forecasted
using either an ARIMA model (\code{method = "arima"}), an AR model (\code{method = "ar"}), 
an exponential smoothing method (\code{method = "ets"}), a linear exponential smoothing 
method allowing missing values (\code{method = "ets.na"}), or a random walk with drift model 
(\code{method = "rwdrift"}). The forecast coefficients are then multiplied by the principal 
components to obtain a forecast curve.}
\usage{
\method{forecast}{ftsm}(object, h = 10, method = c("ets", "arima", "ar", "ets.na", 
 "rwdrift", "rw", "struct", "arfima"), level = 80, jumpchoice = c("fit", 
  "actual"), pimethod = c("parametric", "nonparametric"), B = 100, 
   usedata = nrow(object$coeff), adjust = TRUE, model = NULL,
    damped = NULL, stationary = FALSE, ...)
}
\arguments{
\item{object}{Output from \code{\link[ftsa]{ftsm}}.}
\item{h}{Forecast horizon.}
\item{method}{Univariate time series forecasting methods. Current possibilities are \dQuote{ets}, \dQuote{arima}, \dQuote{ets.na}, \dQuote{rwdrift} and \dQuote{rw}.}
\item{level}{Coverage probability of prediction intervals.}
\item{jumpchoice}{Jump-off point for forecasts. Possibilities are \dQuote{actual} and \dQuote{fit}.
	If \dQuote{actual}, the forecasts are bias-adjusted by the difference between the fit and the last year of observed data.
	Otherwise, no adjustment is used. See Booth et al. (2006) for the detail on jump-off point.}
\item{pimethod}{Indicates if parametric method is used to construct prediction intervals.}
\item{B}{Number of bootstrap samples.}
\item{usedata}{Number of time periods to use in forecasts. Default is to use all.}
\item{adjust}{If \code{adjust = TRUE}, adjusts the variance so that the one-step forecast variance matches the empirical one-step forecast variance.}
\item{model}{If the \code{ets} method is used, \code{model} allows a model specification to be passed to \code{\link[forecast]{ets}()}.}
\item{damped}{If the \code{ets} method is used, \code{damped} allows the damping specification to be passed to \code{\link[forecast]{ets}()}.}
\item{stationary}{If \code{stationary = TRUE}, \code{method} is set to \code{method = "ar"} and only stationary AR models are used.}
\item{...}{Other arguments passed to forecast routine.}
}

\details{
1. Obtain a smooth curve \eqn{f_t(x)} for each \eqn{t} using a nonparametric smoothing technique.

2. Decompose the smooth curves via a functional principal component analysis.

3. Fit a univariate time series model to each of the principal component scores.

4. Forecast the principal component scores using the fitted time series models.

5. Multiply the forecast principal component scores by fixed principal components to obtain forecasts of \eqn{f_{n+h}(x)}.

6. The estimated variances of the error terms (smoothing error and model residual error) are used to compute prediction intervals for the forecasts.
}

\value{List with the following components:
\item{mean}{An object of class \code{fts} containing point forecasts.}
\item{lower}{An object of class \code{fts} containing lower bound for prediction intervals.}
\item{upper}{An object of class \code{fts} containing upper bound for prediction intervals.}
\item{fitted}{An object of class \code{fts} of one-step-ahead forecasts for historical data.}
\item{error}{An object of class \code{fts} of one-step-ahead errors for historical data.}
\item{coeff}{List of objects of type \code{forecast} containing the coefficients and their forecasts.}
\item{coeff.error}{One-step-ahead forecast errors for each of the coefficients.}
\item{var}{List containing the various components of variance: model, error, mean, total and coeff.}
\item{model}{Fitted \code{\link[ftsa]{ftsm}} model.}
\item{bootsamp}{An array of \eqn{dimension = c(p, B, h)} containing the bootstrapped point forecasts. 
                \eqn{p} is the number of variables. \eqn{B} is the number of bootstrap samples.
                \eqn{h} is the forecast horizon.}
}

\references{
H. Booth and R. J. Hyndman and L. Tickle and P. D. Jong (2006) "Lee-Carter mortality forecasting: A multi-country comparison of variants and extensions", \emph{Demographic Research}, \bold{15}, 289-310.

B. Erbas and R. J. Hyndman and D. M. Gertig (2007) "Forecasting age-specific breast cancer mortality using functional data model", \emph{Statistics in Medicine}, \bold{26}(2), 458-470.

R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series" (with discussion), \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.

H. L. Shang (2012) "Functional time series approach for forecasting very short-term electricity demand", \emph{Journal of Applied Statistics}, \bold{40}(1), 152-168.

H. L. Shang (2013) "ftsa: An R package for analyzing functional time series", \emph{The R Journal}, \bold{5}(1), 64-72.

H. L. Shang, A. Wisniowski, J. Bijak, P. W. F. Smith and J. Raymer (2014) "Bayesian functional models for population forecasting", in M. Marsili and G. Capacci (eds), Proceedings of the Sixth Eurostat/UNECE Work Session on Demographic Projections, Istituto nazionale di statistica, Rome, pp. 313-325.

H. L. Shang (2015) "Selection of the optimal Box-Cox transformation parameter for modelling and forecasting age-specific fertility", \emph{Journal of Population Research}, \bold{32}(1), 69-79.

H. L. Shang (2015) "Forecast accuracy comparison of age-specific mortality and life expectancy: Statistical tests of the results", \emph{Population Studies}, \bold{69}(3), 317-335.

H. L. Shang, P. W. F. Smith, J. Bijak, A. Wisniowski (2016) "A multilevel functional data method for forecasting population, with an application to the United Kingdom", \emph{International Journal of Forecasting}, \bold{32}(3), 629-649.
}

\author{Rob J Hyndman}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecastfplsr}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}}
\examples{
# ElNino is an object of class sliced functional time series.
# Via functional principal component decomposition, the dynamic was captured 
# by a few principal components and principal component scores. 
# By using an exponential smoothing method, 
# the principal component scores are forecasted.
# The forecasted curves are constructed by forecasted principal components 
# times fixed principal components plus the mean function.	
forecast(object = ftsm(ElNino_ERSST_region_1and2), h = 10, method = "ets")              
forecast(object = ftsm(ElNino_ERSST_region_1and2, weight = TRUE))
}
\keyword{models}


