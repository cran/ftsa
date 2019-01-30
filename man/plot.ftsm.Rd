\name{plot.ftsm}
\alias{plot.ftsm}
\title{Plot fitted model components for a functional time series model}
\description{
Plot showing the basis functions in the top row of plots and the coefficients in the bottom row of plots. 
}
\usage{
\method{plot}{ftsm}(x, components, components.start = 0, xlab1 = x$y$xname, ylab1 = "Basis function", 
 xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean", 
  level.lab = "Level", main.title = "Main effects", 
   interaction.title = "Interaction", basiscol = 1, coeffcol = 1, 
    outlier.col = 2, outlier.pch = 19, outlier.cex = 0.5, ...)   
}
\arguments{
\item{x}{Output from \code{\link[ftsa]{ftsm}}.}
\item{components}{Number of principal components to plot.}
\item{components.start}{Plotting specified component.}
\item{xlab1}{x-axis label for basis functions.}
\item{xlab2}{x-axis label for coefficient time series.}
\item{ylab1}{y-axis label for basis functions.}
\item{ylab2}{y-axis label for coefficient time series.}
\item{mean.lab}{Label for mean component.}
\item{level.lab}{Label for level component.}
\item{main.title}{Title for main effects.}
\item{interaction.title}{Title for interaction terms.} 
\item{basiscol}{Colors for basis functions if plot.type="components".}
\item{coeffcol}{Colors for time series coefficients if plot.type="components".}
\item{outlier.col}{Colour for outlying years.}
\item{outlier.pch}{Plotting character for outlying years.}
\item{outlier.cex}{Size of plotting character for outlying years.}
\item{...}{Plotting parameters.}
}
\value{
None. Function produces a plot.
}

\references{
R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series" (with discussion), \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{Rob J Hyndman}
\seealso{
\code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}
}
\examples{
# plot different principal components.	
plot.ftsm(ftsm(y = ElNino_ERSST_region_1and2, order = 2), components = 2)
}
\keyword{models}
