\name{FANOVA}
\alias{FANOVA}
\title{
Functional analysis of variance fitted by means. 
}
\description{
Decomposition by functional analysis of variance fitted by means.
}
\usage{
FANOVA(data_pop1, data_pop2, year=1959:2020, age= 0:100, 
	       n_prefectures=51, n_populations=2)
}
\arguments{
  \item{data_pop1}{It's a p by n matrix}
  \item{data_pop2}{It's a p by n matrix}
  \item{year}{Vector with the years considered in each population. }
  \item{n_prefectures}{Number of prefectures}
  \item{age}{Vector with the ages considered in each year.}
  \item{n_populations}{Number of populations.}
}
\value{
\item{FGE_mean}{FGE_mean, a vector of dimension p}
\item{FRE_mean}{FRE_mean, a matrix of dimension length(row_partition_index) by p.}
\item{FCE_mean}{FCE_mean, a matrix of dimension length(column_partition_index) by p.}
}
\references{
C. F. Jimenez Varon, Y. Sun and H. L. Shang (2023) ``Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality".

Ramsay, J. and B. Silverman (2006). Functional Data Analysis. Springer Series in Statistics. Chapter 13. New York: Springer
}
\author{
Cristian Felipe Jimenez Varon, Ying Sun, Han Lin Shang
}

\seealso{
\code{\link{Two_way_median_polish}}
}
\examples{
# The US mortality data  1959-2020 for two populations and three states 
# (New York, California, Illinois)
# Compute the functional Anova decomposition fitted by means.
FANOVA_means <- FANOVA(data_pop1 = t(all_hmd_male_data), 
					      data_pop2 = t(all_hmd_female_data),
					      year = 1959:2020, age =  0:100, 
					      n_prefectures = 3, n_populations = 2)

##1. The funcional grand effect
FGE = FANOVA_means$FGE_mean
##2. The funcional row effect
FRE = FANOVA_means$FRE_mean
##3. The funcional column effect
FCE = FANOVA_means$FCE_mean
}
\keyword{methods}
