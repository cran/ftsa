\name{skew_t_fun}
\alias{skew_t_fun}
\title{
Skewed t distribution
}
\description{
Fitting a parametric skewed t distribution of Fernandez and Steel's (1998) method
}
\usage{
skew_t_fun(data, gridpoints, M = 5001)
}
\arguments{
  \item{data}{a data matrix of dimension \code{n} by \code{p}}
  \item{gridpoints}{Grid points}
  \item{M}{number of grid points}
}
\details{
1) Fit a skewed t distribution to data, and obtain four latent parameters;
2) Transform the four latent parameters so that they are un-constrained;
3) Fit a vector autoregressive model to these transformed latent parameters;
4) Obtain their forecasts, and then back-transform them to the original scales;
5) Via the skewed t distribution in Step 1), we obtain forecast density using the forecast latent parameters.
}
\value{
\item{m}{Grid points within data range}
\item{skewed_t_den_fore}{Density forecasts via a skewed t distribution}
}
\references{
Fernandez, C. and Steel, M. F. J. (1998), `On Bayesian modeling of fat tails and skewness', \emph{Journal of the American Statistical Association: Theory and Methods}, \bold{93}(441), 359-371.
}
\author{
Han Lin Shang
}
\note{
This is a parametric approach for fitting and forecasting density.
}
\seealso{
\code{\link[ftsa]{CoDa_FPCA}}, \code{\link[ftsa]{Horta_Ziegelmann_FPCA}}, \code{\link[ftsa]{LQDT_FPCA}}
}
\examples{
skew_t_fun(DJI_return)
}
\keyword{distribution}
\keyword{models}
