\name{Two_way_Residuals_means}
\alias{Two_way_Residuals_means}
\title{
Functional time series decomposition into deterministic (functional analysis of variance fitted by means), and time-varying components (functional residuals).
}
\description{
Decomposition of functional time series into deterministic (by functional analysis of variance fitted by means), and time-varying components (functional residuals)
}
\usage{
Two_way_Residuals_means(data_pop1, data_pop2, year, age, n_prefectures, n_populations)
}
\arguments{
  \item{data_pop1}{A p by n matrix}
  \item{data_pop2}{A p by n matrix}
  \item{year}{Vector with the years considered in each population. }
  \item{n_prefectures}{Number of prefectures}
  \item{age}{Vector with the ages considered in each year.}
  \item{n_populations}{Number of populations.}
}
\value{
\item{residuals1}{A matrix with  dimension n by p.}
\item{residuals2}{A matrix with  dimension n by p.}
\item{rd}{A two dimension logic vector proving that the decomposition sum up the data.}
\item{R}{A matrix of dimension as n by 2p. This represents the time-varying component in the decomposition.}
\item{Fixed_comp}{A matrix of dimension as n by 2p. This represents the deterministic component in the decomposition.}
}
\references{
C. F. Jimenez Varon, Y. Sun and H. L. Shang (2023) ``Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality".

Ramsay, J. and B. Silverman (2006). Functional Data Analysis. Springer Series in Statistics. Chapter 13. New York: Springer.
}
\author{
Cristian Felipe Jimenez Varon, Ying Sun, Han Lin Shang
}
\seealso{
\code{\link{Two_way_Residuals}}
}
\examples{
# The US mortality data  1959-2020, for two populations
# and three states (New York, California, Illinois)
# Compute the functional Anova decomposition fitted by means.
FANOVA_means_residuals <- Two_way_Residuals_means(data_pop1=t(all_hmd_male_data),
                            data_pop2=t(all_hmd_female_data), year = 1959:2020,
                            age = 0:100, n_prefectures = 3, n_populations = 2)
                            
# The results
##1. The functional residuals from population 1
Residuals_pop_1=FANOVA_means_residuals$residuals1
##2. The functional residuals from population 2
Residuals_pop_2=FANOVA_means_residuals$residuals2
##3. A logic vector whose components indicate whether the sum of deterministic
##  and time-varying components recover the original FTS.
Construct_data=FANOVA_means_residuals$rd
##4. Time-varying components for all the populations. The functional residuals
All_pop_functional_residuals <- FANOVA_means_residuals$R
##5. The deterministic components from the functional ANOVA decomposition
deterministic_comp <- FANOVA_means_residuals$Fixed_comp
}
\keyword{methods}
