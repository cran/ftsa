\name{CoDa_BayesNW}
\alias{CoDa_BayesNW}
\title{
Compositional data analytic approach and nonparametric function-on-function regression for forecasting density
}
\description{
Log-ratio transformation from constrained space to unconstrained space, where a standard nonparametric function-on-function regression can be applied.
}
\usage{
CoDa_BayesNW(data, normalization, m = 5001, 
	band_choice = c("Silverman", "DPI"), 
	kernel = c("gaussian", "epanechnikov"))
}
\arguments{
  \item{data}{Densities or raw data matrix of dimension N by p, where N denotes sample size and p denotes dimensionality}
  \item{normalization}{If a standardization should be performed?}
  \item{m}{Grid points within the data range}
  \item{band_choice}{Selection of optimal bandwidth}
  \item{kernel}{Type of kernel function}
}
\details{
1) Compute the geometric mean function
2) Apply the centered log-ratio transformation
3) Apply a nonparametric function-on-function regression to the transformed data
4) Transform forecasts back to the compositional data
5) Add back the geometric means, to obtain the forecasts of the density function
}
\value{
Out-of-sample density forecasts
}
\references{
Egozcue, J. J., Diaz-Barrero, J. L. and Pawlowsky-Glahn, V. (2006) `Hilbert space of probability density functions based on Aitchison geometry', \emph{Acta Mathematica Sinica}, \bold{22}, 1175-1182.

Ferraty, F. and Shang, H. L. (2021) `Nonparametric density-on-density regression', working paper.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{CoDa_FPCA}}
}
\examples{
\dontrun{
CoDa_BayesNW(data = DJI_return, normalization = "TRUE", 
		band_choice = "DPI", kernel = "epanechnikov")
}	
}
\keyword{methods}
