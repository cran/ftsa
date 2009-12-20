\name{dynupdateweightselect}
\Rdversion{1.1}
\alias{dynupdateweightselect}
\title{
Selecting an optimal tuning parameter for the PLS and RR methods
}
\description{
Selecting an optimal tuning parameter for the penalized least squares and ridge regression methods used in the \code{\link[ftsa]{dynupdate}}. 
}
\usage{
dynupdateweightselect(data, method = c("pls", "ridge"), interval = c(0, 10^4), p, 
                      backh = 1, errortype = c("mse", "mae", "mape"))
}
\arguments{
  \item{data}{An object of class \code{sfts}.}
  \item{method}{Either the penalized least squares method or the ridge regression method.}
  \item{interval}{Interval where the optimization takes place.}
  \item{p}{The most recent time period within the function support range.}
  \item{backh}{The number of curves used as the validation data. For instance, when \code{backh = 2}, the last two curves are subsequently treated as the validation data, while the remaining curves are treated as the training data.}
  \item{errortype}{Type of error measure.}
}
\details{
The optimal tuning parameter is selected by minimizing an averaged forecast error measure. This function also utilizes the \code{\link[stats]{optimize}} function. 
}
\value{
\item{minimum}{Optimal tuning parameter.}
\item{objective}{Value of the error function at the optimal point.}
}
\references{
H. Shang and R. J. Hyndman (2009) "Nonparametric time series forecasting with dynamic updating", In R. S. Anderssen, R. D. Braddock and L.T.H. Newham (eds), 18th World IMACS Congress and MODSIM09 International Congress on Modelling and Simulation. Modelling and Simulation Society of Australia and New Zealand and International Association for Mathematics and Computers in Simulation, July 2009, pp. 1552-1558. ISBN: 978-0-9758400-7-8. 
\url{http://www.mssanz.org.au/modsim09/D11/shang.pdf}
}
\author{
Han Lin Shang
}
\note{
This function can be computationally intensive for \code{backh > 1}. 
}
\seealso{
\code{\link[ftsa]{ftsmweightselect}}, \code{\link[ftsa]{fplsrweightselect}}, \code{\link[ftsa]{dynupdate}}
}
\examples{
dynupdateweightselect(data = ElNino, method = "pls", p = 4, errortype = "mse")
}
\keyword{models}

