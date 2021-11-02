\name{CoDa_FPCA}
\alias{CoDa_FPCA}
\title{
Compositional data analytic approach and functional principal component analysis for forecasting density
}
\description{
Log-ratio transformation from constrained space to unconstrained space, where a standard functional principal component analysis can be applied.
}
\usage{
CoDa_FPCA(data, normalization, h_scale = 1, m = 5001, 
	band_choice = c("Silverman", "DPI"), 
	kernel = c("gaussian", "epanechnikov"), 
	varprop = 0.99, fmethod)
}
\arguments{
  \item{data}{Densities or raw data matrix of dimension n by p, where n denotes sample size and p denotes dimensionality}
  \item{normalization}{If a standardization should be performed?}
  \item{h_scale}{Scaling parameter in the kernel density estimator}
  \item{m}{Grid point within the data range}
  \item{band_choice}{Selection of optimal bandwidth}
  \item{kernel}{Type of kernel functions}
  \item{varprop}{Proportion of variance explained}
  \item{fmethod}{Univariate time series forecasting method}
}
\details{
1) Compute the geometric mean function
2) Apply the centered log-ratio transformation
3) Apply FPCA to the transformed data
4) Forecast principal component scores
5) Transform forecasts back to the compositional data
6) Add back the geometric means, to obtain the forecasts of the density function
}
\value{
Out-of-sample forecast densities
}
\references{
Boucher, M.-P. B., Canudas-Romo, V., Oeppen, J. and Vaupel, J. W. (2017) `Coherent forecasts of mortality with compositional data analysis', \emph{Demographic Research}, \bold{37}, 527-566.

Egozcue, J. J., Diaz-Barrero, J. L. and Pawlowsky-Glahn, V. (2006) `Hilbert space of probability density functions based on Aitchison geometry', \emph{Acta Mathematica Sinica}, \bold{22}, 1175-1182.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{Horta_Ziegelmann_FPCA}}, \code{\link[ftsa]{LQDT_FPCA}}, \code{\link[ftsa]{skew_t_fun}}
}
\examples{
\dontrun{
CoDa_FPCA(data = DJI_return, normalization = "TRUE", band_choice = "DPI", 
	kernel = "epanechnikov", varprop = 0.9, fmethod = "ETS")
}
}
\keyword{methods}
