\name{Two_way_Residuals}
\alias{Two_way_Residuals}
\title{
Functional time series decomposition into deterministic (from functional median polish from Sun and Genton (2012)), and time-varying components (functional residuals).
}
\description{
Decomposition of functional time series into deterministic (from functional median polish), and time-varying components (functional residuals)
}
\usage{
Two_way_Residuals(Y, n_prefectures, year, age, n_populations)
}
\arguments{
  \item{Y}{A matrix with dimension n by 2p. The functional data}
  \item{year}{Vector with the years considered in each population}
  \item{n_prefectures}{Number of prefectures}
  \item{age}{Vector with the ages considered in each year}
  \item{n_populations}{Number of populations}
}
\value{
\item{residuals1 }{A matrix with dimension n by p}
\item{residuals2 }{A matrix with dimension n by p}
\item{rd}{A two dimension logic vector that proves that the decomposition sum up to the data}
\item{R}{A matrix with the same dimension as Y. This represent the time-varying component in the decomposition}
\item{Fixed_comp}{A matrix with the same dimension as Y. This represent the deterministic component in the decomposition}
}
\references{
C. F. Jimenez Varon, Y. Sun and H. L. Shang (2023) "Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality".

Sun, Ying, and Marc G. Genton (2012). "Functional Median Polish". Journal of Agricultural, Biological, and Environmental Statistics 17(3), 354-376.
}
\author{
Cristian Felipe Jimenez Varon, Ying Sun, Han Lin Shang
}
\seealso{
\code{\link{Two_way_Residuals_means}}
}
\examples{
# The US mortality data  1959-2020, for two populations
# and three states (New York, California, Illinois)
# Column binds the data from both populations
Y = cbind(all_hmd_male_data, all_hmd_female_data)
# Decompose FTS into deterministic (from functional median polish)
# and time-varying components (functional residuals).
FMP_residuals <- Two_way_Residuals(Y,n_prefectures=3,year=1959:2020,
                                   age=0:100,n_populations=2)
# The results
##1. The functional residuals from population 1
Residuals_pop_1=FMP_residuals$residuals1
##2. The functional residuals from population 2
Residuals_pop_2=FMP_residuals$residuals2
##3. A logic vector whose components indicate whether the sum of deterministic
##   and time-varying components recover the original FTS.
Construct_data=FMP_residuals$rd
##4. Time-varying components for all the populations. The functional residuals
All_pop_functional_residuals <- FMP_residuals$R
##5. The deterministic components from the functional median polish decomposition
deterministic_comp <- FMP_residuals$Fixed_comp
}
\keyword{methods}
