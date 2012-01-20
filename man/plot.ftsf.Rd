\name{plot.ftsf}
\alias{plot.ftsf}
\title{Plot fitted model components for a functional time series model}
\description{
Plot fitted model components for a \code{fts} object.
}
\usage{
\method{plot}{ftsf}(x, plot.type = c("function", "components", "variance"), 
 components, xlab1 = fit$y$xname, ylab1 = "Basis function", 
  xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean", 
   level.lab = "Level", main.title = "Main effects", 
    interaction.title = "Interaction", vcol = 1:3, shadecols = 7, 
     fcol = 4, basiscol = 1, coeffcol = 1, outlier.col = 2,
      outlier.pch = 19, outlier.cex = 0.5,...)
}
\arguments{
\item{x}{Output from \code{\link[ftsa]{forecast.ftsm}}.}
\item{plot.type}{Type of plot.}
\item{components}{Number of principal components.}
\item{xlab1}{x-axis label for principal components.}
\item{xlab2}{x-axis label for coefficient time series.}
\item{ylab1}{y-axis label for principal components.}
\item{ylab2}{y-axis label for coefficient time series.}
\item{mean.lab}{Label for mean component.}
\item{level.lab}{Label for level component.}
\item{main.title}{Title for main effects.}
\item{interaction.title}{Title for interaction terms.}
\item{vcol}{Colors to use if \code{plot.type = "variance"}.}
\item{shadecols}{Color for shading of prediction intervals when \code{plot.type = "components"}.}
\item{fcol}{Color of point forecasts when \code{plot.type = "components"}.}
\item{basiscol}{Colors for principal components if \code{plot.type = "components"}.}
\item{coeffcol}{Colors for time series coefficients if \code{plot.type = "components"}.}
\item{outlier.col}{Colors for outlying years.}
\item{outlier.pch}{Plotting character for outlying years.}
\item{outlier.cex}{Size of plotting character for outlying years.}
\item{...}{Plotting parameters.} 
}
\details{
When \code{plot.type = "function"}, it produces a plot of the forecast functions; 

When \code{plot.type = "components"}, it produces a plot of the principla components and coefficients with forecasts and prediction intervals for each coefficient;

When \code{plot.type = "variance"}, it produces a plot of the variance components.
}
\value{Function produces a plot.}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.fmres}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}}
\author{Rob J Hyndman}
\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\examples{
plot(x = forecast(object = ftsm(y = ElNino)))
}
\keyword{hplot}

