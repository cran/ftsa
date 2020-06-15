\name{isfe.fts}
\alias{isfe.fts}
\title{Integrated Squared Forecast Error for models of various orders}
\description{
Computes integrated squared forecast error (ISFE) values for functional time series models of various orders.
}
\usage{
isfe.fts(data, max.order = N - 3, N = 10, h = 5:10, method = 
 c("classical", "M", "rapca"), mean = TRUE, level = FALSE, 
  fmethod = c("arima", "ar", "ets", "ets.na", "struct", "rwdrift", 
   "rw", "arfima"), lambda = 3, ...)
}
\arguments{
\item{data}{An object of class \code{fts}.}
\item{max.order}{Maximum number of principal components to fit.}
\item{N}{Minimum number of functional observations to be used in fitting a model.}
\item{h}{Forecast horizons over which to average.}
\item{method}{Method to use for principal components decomposition. Possibilities are \dQuote{M}, \dQuote{rapca} and \dQuote{classical}.}
\item{mean}{Indicates if mean term should be included.}
\item{level}{Indicates if level term should be included.}
\item{fmethod}{Method used for forecasting. Current possibilities are \dQuote{ets}, \dQuote{arima}, \dQuote{ets.na},
	\dQuote{struct}, \dQuote{rwdrift} and \dQuote{rw}.}
\item{lambda}{Tuning parameter for robustness when \code{method = "M"}.}
\item{...}{Additional arguments controlling the fitting procedure.}
}
\value{
Numeric matrix with \code{(max.order+1)} rows and \code{length(h)} columns
containing ISFE values for models of orders \code{0:(max.order)}. 
}
\author{Rob J Hyndman}
\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.
}
\note{
This function can be very time consuming for data with large dimensionality or large sample size.
By setting \code{max.order} small, computational speed can be dramatically increased.
}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.fmres}}, \code{\link[ftsa]{summary.fm}}, \code{\link[ftsa]{residuals.fm}}}
\keyword{models}

