\name{forecastfplsr}
\alias{forecastfplsr}
\title{Forecast functional time series}
\description{
The decentralized response is forecasted by multiplying the estimated regression coefficient with the new decentralized predictor
}
\usage{
forecastfplsr(object, components, h)
}
\arguments{
  \item{object}{An object of class \code{fts}.}
  \item{components}{Number of optimal components.}
  \item{h}{Forecast horizon.}
}
\value{
A \code{fts} class object, containing forecasts of responses.
}
\references{
R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series" (with discussion), \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}
}
\examples{
# A set of functions are decomposed by functional partial least squares decomposition.	
# By forecasting univariate partial least squares scores, the forecasted curves are 
# obtained by multiplying the forecasted scores by fixed functional partial least 
# squares function plus fixed mean function.
forecastfplsr(object = ElNino_ERSST_region_1and2, components = 2, h = 5)
}
\keyword{models}

