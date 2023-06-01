\name{Two_way_median_polish}
\alias{Two_way_median_polish}
\title{
Two-way functional median polish from Sun and Genton (2012)
}
\description{
Decomposition by two-way functional median polish
}
\usage{
Two_way_median_polish(Y, year=1959:2020, age=0:100, n_prefectures=51, n_populations=2)
}
\arguments{
  \item{Y}{A  matrix with dimension n by 2p. The functional data.}
  \item{year}{Vector with the years considered in each population. }
  \item{n_prefectures}{Number of prefectures}
  \item{age}{Vector with the ages considered in each year.}
  \item{n_populations}{Number of populations.}
}
\value{
\item{grand_effect}{grand_effect, a vector of dimension p}
\item{row_effect}{row_effect, a matrix of dimension length(row_partition_index) by p.}
\item{col_effect}{col_effect, a matrix of dimension length(column_partition_index) by p}
}
\references{
C. F. Jim\'{e}nez Var\'{o}n, Y. Sun and H. L. Shang (2023) ``Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality".

Sun, Ying, and Marc G. Genton (2012) ``Functional Median Polish", Journal of Agricultural, Biological, and Environmental Statistics, 17(3), 354-376. 
}

\author{
Cristian Felipe Jim\'{e}nez Var\'{o}n, Ying Sun, Han Lin Shang
}
\seealso{
\code{\link{FANOVA}}
}
\examples{
# The US mortality data  1959-2020 for two populations and three states 
# (New York, California, Illinois)
# Compute the functional median polish decomposition.
FMP = Two_way_median_polish(cbind(all_hmd_male_data, all_hmd_female_data), 
		n_prefectures = 3, year = 1959:2020, age = 0:100, n_populations = 2)

##1. The funcional grand effect
FGE = FMP$grand_effect
##2. The funcional row effect
FRE = FMP$row_effect
##3. The funcional column effect
FCE = FMP$col_effect
}
\keyword{methods}
