\name{plotfplsr}
\alias{plotfplsr}
\title{
Plot fitted model components for a functional time series model
}
\description{
Plot showing the basis functions of the predictors in the top row, followed by the basis functions of the responses in the second row, then the coefficients in the bottom row of plots. 
}
\usage{
plotfplsr(x, xlab1 = x$ypred$xname, ylab1 = "Basis function", xlab2 = "Time", 
 ylab2 = "Coefficient", mean.lab = "Mean", interaction.title = "Interaction")
}
\arguments{
  \item{x}{Output from \code{\link[ftsa]{fplsr}}.}
  \item{xlab1}{x-axis label for basis functions.}
  \item{ylab1}{y-axis label for basis functions.}
  \item{xlab2}{x-axis label for coefficient time series.}
  \item{ylab2}{y-axis label for coefficient time series.}
  \item{mean.lab}{Label for mean component.}
  \item{interaction.title}{Title for interaction terms.}
}
\value{
None. Function produces a plot.
}
\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series" (with discussion), \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}
}
\examples{
# Fit the data by the functional partial least squares.	
ausfplsr = fplsr(data = ElNino_ERSST_region_1and2, order = 2)
plotfplsr(x = ausfplsr)
}
\keyword{models}

