\name{ftsmweightselect}
\alias{ftsmweightselect}
\title{
Selection of the weight parameter used in the weighted functional time series model.
}
\description{
The geometrically decaying weights are used to estimate the mean curve and functional principal components, where more weights are assigned to the more recent data than the data from the distant past.
}
\usage{
ftsmweightselect(data, ncomp = 6, ntestyear, errorcriterion = c("mae", "mse", "mape"))
}
\arguments{
  \item{data}{An object of class \code{fts}.}
  \item{ncomp}{Number of components.}
  \item{ntestyear}{Number of holdout observations used to assess the forecast accuracy.}
  \item{errorcriterion}{Error measure.}
}
\details{
The data set is split into a fitting period and forecasting period. Using the data in the fitting period, we compute the one-step-ahead forecasts and calculate
the forecast error. Then, we increase the fitting period by one, and carry out the same forecasting procedure until the fitting period covers entire data set. 
The forecast accuracy is determined by the averaged forecast error across the years in the forecasting period. By using an optimization algorithm, we select the
optimal weight parameter that would result in the minimum forecast error.
}
\value{
Optimal weight parameter.
}
\references{
R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{
Han Lin Shang
}
\note{
Can be computational intensive, as it takes about half-minute to compute. For example, ftsmweightselect(ElNinosmooth, ntestyear = 1).
}
\seealso{
\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}
}
\keyword{methods}

