\name{plot.fm}
\alias{plot.fm}
\title{Plot fitted model components for a functional model}
\description{
When \code{class(x)[1] = ftsm}, plot showing the principal components in the top row of plots and the coefficients in the bottom row of plots.

When \code{class(x)[1] = fm}, plot showing the predictor scores in the top row of plots and the response loadings in the bottom row of plots.
}
\usage{
\method{plot}{fm}(x, order, xlab1 = x$y$xname, ylab1 = "Principal component", 
 xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean", 
  level.lab = "Level", main.title = "Main effects", interaction.title 
   = "Interaction", basiscol = 1, coeffcol = 1, outlier.col = 2, 
    outlier.pch = 19, outlier.cex = 0.5, ...)
}
\arguments{
\item{x}{Output from \code{\link[ftsa]{ftsm}} or \code{\link[ftsa]{fplsr}}.}
\item{order}{Number of principal components to plot. Default is all principal components in a model.}
\item{xlab1}{x-axis label for principal components.}
\item{xlab2}{x-axis label for coefficient time series.}
\item{ylab1}{y-axis label for principal components.}
\item{ylab2}{y-axis label for coefficient time series.}
\item{mean.lab}{Label for mean component.}
\item{level.lab}{Label for level component.}
\item{main.title}{Title for main effects.}
\item{interaction.title}{Title for interaction terms.}
\item{basiscol}{Colors for principal components if \code{plot.type = "components"}.}
\item{coeffcol}{Colors for time series coefficients if \code{plot.type = "components"}.}
\item{outlier.col}{Colors for outlying years.}
\item{outlier.pch}{Plotting character for outlying years.}
\item{outlier.cex}{Size of plotting character for outlying years.}
\item{...}{Plotting parameters.}
}
\value{
Function produces a plot.
}
\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}, \code{\link[ftsa]{plot.fmres}}, \code{\link[ftsa]{plot.ftsf}}}
\author{Rob J Hyndman}
\examples{
plot(x = ftsm(y = ElNino_ERSST_region_1and2))
}
\keyword{hplot}

