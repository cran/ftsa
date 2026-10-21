\name{pm_10_GR}
\alias{pm_10_GR}
\alias{pm_10_GR_sqrt}
\docType{data}
\title{
Particulate Matter Concentrations (pm10)
}
\description{
This data set consists of half-hourly measurements of the concentrations (measured in ug/m3) of particulate matter with
an aerodynamic diameter of less than 10 \eqn{\mu}m, abbreviated PM10, in ambient air taken in Graz-Mitte, Austria from October 1, 
2010 until March 31, 2011. To stabilize the variance, a square-root transformation can be applied to the data.
}
\usage{
data(pm_10_GR)
}
\format{
A data set.
}
\details{
As epidemiological and toxicological studies have pointed to negative health effects, European Union (EU) regulation sets
pollution standards for the level of the concentration. Policy makers have to ensure compliance with these EU rules and need 
reliable statistical tools to determine appropriate measures, such as partial traffic regulation, and to justify them to the public (see Stadlober, H\"ormann and Pfeiler, 2008). 
}
\source{
Thanks to Professor Siegfried H\"ormann for providing this data set. The original data source is \url{https://zenodo.org/records/7959116} 
}
\references{
A. Aue, D. D. Norinho, S. H\"ormann (2015) "On the prediction of stationary functional time series", \emph{Journal of the American Statistical Association}, \bold{110}(509), 378-392.

E. Stadlober, S. H\"ormann, B. Pfeiler (2008) "Quality and performance of a PM10 daily forecasting model", \emph{Atmospheric Environment}, \bold{42}, 1098-1109.

S. H\"ormann, B. Pfeiler, E. Stadlober (2005) "Analysis and prediction of particulate matter PM10 for the winter season in Graz", \emph{Austrian Journal of Statistics}, \bold{34}(4), 307-326.

H. L. Shang (2017) "Functional time series forecasting with dynamic updating: An application to intraday particulate matter concentration", \emph{Econometrics and Statistics}, \bold{1}, 184-200. 
}
\examples{
plot(pm_10_GR)
}
\keyword{datasets}
