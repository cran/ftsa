\name{residuals.fm}
\alias{residuals.fm}
\title{Compute residuals from a functional model}
\description{After fitting a functional model, it is useful to inspect the residuals. 
This function extracts the relevant information from the fit object and puts it in a form suitable for
plotting with \code{image}, \code{persp}, \code{contour}, \code{filled.contour}, etc.}

\usage{
\method{residuals}{fm}(object, ...)
}
\arguments{
\item{object}{Output from \code{\link[ftsa]{ftsm}} or \code{\link[ftsa]{fplsr}}.}
\item{...}{Other arguments.}
}
\value{
Produces an object of class \dQuote{fmres} containing the residuals from the model.
}
\references{
B. Erbas and R. J. Hyndman and D. M. Gertig (2007) "Forecasting age-specific breast cancer mortality using functional data model", \emph{Statistics in Medicine}, \bold{26}(2), 458-470.

R. J. Hyndman and M. S. Ullah (2007) "Robust forecasting of mortality and fertility rates: A functional data approach", \emph{Computational Statistics \& Data Analysis}, \bold{51}(10), 4942-4956.

R. J. Hyndman and H. Booth (2008) "Stochastic population forecasts using functional data models for mortality, fertility and migration", \emph{International Journal of Forecasting}, \bold{24}(3), 323-342.
}
\author{Rob J Hyndman}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{summary.fm}}, \code{\link[ftsa]{plot.fm}}, \code{\link[ftsa]{plot.fmres}}}
\examples{
plot(residuals(object = ftsm(ElNino)), xlab = "Year", ylab = "Month")
}
\keyword{models}
