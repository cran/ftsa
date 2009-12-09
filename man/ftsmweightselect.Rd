\name{ftsmweightselect}
\Rdversion{1.1}
\alias{ftsmweightselect}
\title{
Selecting an optimal weight parameter for the weighted functional principal component regression
}
\description{
Selecting an optimal geometrically decaying weight used in the \code{\link[ftsa]{ftsm}}. 
}
\usage{
ftsmweightselect(data, errortype = c("L1norm", "L2norm", "Linfnorm"), h)
}
\arguments{
  \item{data}{An object of class \code{fts}.}
  \item{errortype}{Type of error measure.}
  \item{h}{Forecast horizon.}
}
\details{
The optimal weight parameter is selected by minimizing an error measure, using the \code{\link[stats]{optimize}} function.
}
\value{
\item{minimum}{Optimal weight parameter.}
\item{objective}{Optimal error.}
}
\references{
R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{fplsrweightselect}}
}
\examples{
ftsmweightselect(Australiasmoothfertility, errortype = "L2norm", h = 10)
}
\keyword{models}

