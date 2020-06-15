\name{ftsm}
\alias{ftsm}
\title{Fit functional time series model}
\description{
Fits a principal component model to a \code{fts} object. The
function uses optimal orthonormal principal components obtained from a
principal components decomposition.
}
\usage{
ftsm(y, order = 6, ngrid = max(500, ncol(y$y)), method = c("classical", 
 "M", "rapca"), mean = TRUE, level = FALSE, lambda = 3, 
  weight = FALSE, beta = 0.1, ...)
}
\arguments{
\item{y}{An object of class \code{fts}.}
\item{order}{Number of principal components to fit.}
\item{ngrid}{Number of grid points to use in calculations. Set to maximum of 500 and \code{ncol(y$y)}.}
\item{method}{Method to use for principal components decomposition. Possibilities are \dQuote{M}, \dQuote{rapca} and \dQuote{classical}.}
\item{mean}{If \code{mean = TRUE}, it will estimate mean term in the model before computing basis terms. 
            If \code{mean = FALSE}, the mean term is assumed to be zero.}
\item{level}{If \code{mean = TRUE}, it will include an additional (intercept) term that depends on \eqn{t} but not on \eqn{x}.}
\item{lambda}{Tuning parameter for robustness when \code{method = "M"}.}
\item{weight}{When \code{weight = TRUE}, a set of geometrically decaying weights is applied to the decentralized data.}
\item{beta}{When \code{weight = TRUE}, the speed of geometric decay is governed by a weight parameter.}
\item{...}{Additional arguments controlling the fitting procedure.}
}

\details{
If \code{method = "classical"}, then standard functional principal component decomposition is used, as described by
Ramsay and Dalzell (1991). 

If \code{method = "rapca"}, then the robust principal component algorithm of Hubert, Rousseeuw and Verboven (2002) is used. 

If \code{method = "M"}, then the hybrid algorithm of Hyndman and Ullah (2005) is used.
}

\value{Object of class \dQuote{ftsm} with the following components:
\item{x1}{Time period of a \code{fts} object, which can be obtained from \code{colnames(y$y)}.}
\item{y1}{Variables of a \code{fts} object, which can be obtained from \code{y$x}.}
\item{y}{Original functional time series or sliced functional time series.}
\item{basis}{Matrix of principal components evaluated at value of \code{y$x} (one column for each principal component).
    The first column is the fitted mean or median.}
\item{basis2}{Matrix of principal components excluded from the selected model.}
\item{coeff}{Matrix of coefficients (one column for each coefficient series). The first column is all ones.}
\item{coeff2}{Matrix of coefficients associated with the principal components excluded from the selected model.}
\item{fitted}{An object of class \code{fts} containing the fitted values.}
\item{residuals}{An object of class \code{fts} containing the regression residuals (difference between observed and fitted).}
\item{varprop}{Proportion of variation explained by each principal component.}
\item{wt}{Weight associated with each time period.}
\item{v}{Measure of variation for each time period.}
\item{mean.se}{Measure of standar error associated with the mean.}
}

\author{Rob J Hyndman}
\references{
J. O. Ramsay and C. J. Dalzell (1991) "Some tools for functional data analysis (with discussion)", \emph{Journal of the Royal Statistical Society: Series B}, \bold{53}(3), 539-572.

M. Hubert and P. J. Rousseeuw and S. Verboven (2002) "A fast robust method for principal components with applications to chemometrics", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{60}(1-2), 101-111.

B. Erbas and R. J. Hyndman and D. M. Gertig (2007) "Forecasting age-specific breast cancer mortality using functional data model", \emph{Statistics in Medicine}, \bold{26}(2), 458-470.

R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics and Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}

\seealso{
\code{\link[ftsa]{ftsmweightselect}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.ftsf}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{summary.fm}}
}
\examples{
# ElNino is an object of class sliced functional time series, constructed 
# from a univariate time series. 
# By default, all observations are assigned with equal weighting. 	
ftsm(y = ElNino_ERSST_region_1and2, order = 6, method = "classical", weight = FALSE)
# When weight = TRUE, geometrically decaying weights are used.
ftsm(y = ElNino_ERSST_region_1and2, order = 6, method = "classical", weight = TRUE)
}
\keyword{models}
