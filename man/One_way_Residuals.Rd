\name{One_way_Residuals}
\alias{One_way_Residuals}
\title{
Functional time series decomposition into deterministic (from functional median polish of Sun and Genton (2012)), and functional residual components.
}
\description{
Decomposition of functional time series into deterministic (from functional median polish), and functional residuals
}
\usage{
One_way_Residuals(Y, n_prefectures = 51, year = 1959:2020, age = 0:100)
}
\arguments{
  \item{Y}{The multivariate functional data, which are a matrix with dimension n by 2p, where n is the sample size and p is the dimensionality.}
  \item{n_prefectures}{Number of prefectures.}
  \item{year}{Vector with the years considered in each population.}
  \item{age}{Vector with the ages considered in each year.}
}
\value{
A matrix of dimension n by p.
}
\references{
C. F. Jimenez Varon, Y. Sun and H. L. Shang (2023) ``Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality", arXiv.
\\
Y. Sun and M. G. Genton (2012) ``Functional median polish", Journal of Agricultural, Biological, and Environmental Statistics, 17(3), 354-376.
}
\author{
Cristian Felipe Jimenez Varon, Ying Sun, Han Lin Shang
}
\seealso{
\code{\link{One_way_median_polish}}
}
\examples{
# The US mortality data  1959-2020, for one populations (female) 
# and 3 states (New York, California, Illinois)
# first define the parameters and the row  partitions.
# Define some parameters.
year = 1959:2020
age = 0:100
n_prefectures = 3

#Load the US data. Make sure it is a matrix. 
Y <- all_hmd_female_data
# The results
# Compute the functional residuals. 
FMP_residuals <- One_way_Residuals(Y, n_prefectures=3, year=1959:2020, age=0:100)
}
\keyword{methods}
