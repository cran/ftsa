\name{all_hmd_female_data}
\alias{all_hmd_female_data}
\docType{data}
\title{
The US female log-mortality rate from 1959-2020 and 3 states (New York, California, Illinois).
}
\description{
We generate for the female population in the US. The functional time series corresponding to the log mortality data in each of the 3 states. Each functional time series comprises the ages from 0 to 100+.
}
\usage{data("all_hmd_male_data")}
\format{
   A n x p matrix with n=186 observations on the following p=101 ages from 0 to 100+.
}
\details{
The data generated corresponds to the FTS for the female US log-mortality. The matrix contains 186 FTS stacked by rows. They correspond to 62 (number of years) times 3 (states). Each FTS contains 101 functional values. 
}
\references{
United States Mortality Database (2023). University of California, Berkeley (USA). Department of
Demography at the University of California, Berkeley. Available at usa.mortality.org (data
downloaded on March 15, 2023).
}
\examples{
data(all_hmd_male_data)
}
\keyword{datasets}
