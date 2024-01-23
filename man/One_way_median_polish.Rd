\name{One_way_median_polish}
\alias{One_way_median_polish}
\title{
One-way functional median polish from Sun and Genton (2012)
}
\description{
Decomposition by one-way functional median polish.
}
\usage{
One_way_median_polish(Y, n_prefectures=51, year=1959:2020, age=0:100)
}
\arguments{
  \item{Y}{The multivariate functional data, which are a  matrix with dimension n by 2p, where n is the sample size and p is the dimensionality.}
 \item{year}{Vector with the years considered in each population.}
  \item{n_prefectures}{Number of prefectures.}
  \item{age}{Vector with the ages considered in each year.}
}
\value{
\item{grand_effect}{Grand_effect, a vector of dimension p.}
\item{row_effect}{Row_effect, a matrix of dimension length(row_partition_index) by p.}
}
\references{
C. F. Jimenez Varon, Y. Sun and H. L. Shang (2023) ``Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality", arXiv.
\\
Sun, Ying, and Marc G. Genton (2012) ``Functional Median Polish", Journal of Agricultural, Biological, and Environmental Statistics 17(3), 354-376. 
}

\author{
Cristian Felipe Jimenez Varon, Ying Sun, Han Lin Shang
}
\seealso{
\code{\link{One_way_Residuals}}, \code{\link{Two_way_median_polish}}, \code{\link{Two_way_Residuals}}
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
Y <-  all_hmd_female_data
# Compute the functional median polish decomposition. 
FMP <- One_way_median_polish(Y,n_prefectures=3,year=1959:2020,age=0:100)
# The results
##1. The funcional grand effect
FGE <- FMP$grand_effect
##2. The funcional row effect
FRE <- FMP$row_effect
}
\keyword{methods}
