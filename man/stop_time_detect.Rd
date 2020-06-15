\name{stop_time_detect}
\alias{stop_time_detect}
\title{
Detection of the optimal stopping time in a curve time series
}
\description{
Detecting the optimal stopping time for the glue curing of wood panels in an automatic process environment. 
}
\usage{
stop_time_detect(data, forecasting_method = c("ets", "arima", "rw"))
}
\arguments{
  \item{data}{An object of class \code{fts}}
  \item{forecasting_method}{A univariate time series forecasting method}
}
\value{
\item{break_points_strucchange}{Breakpoints detected by the regression approach}
\item{break_points_ecp}{Breakpoints detected by the distance-based approach}
\item{err_forward}{Forward integrated squared forecast errors}
\item{err_backward}{Backward integrated squared forecast errors (ISFEs)}
\item{ncomp_select_forward}{Number of components selected by the eigenvalue ratio tests based on the forward ISFEs}
\item{ncomp_select_backward}{Number of components selected by the eigenvalue ratio tests based on the backward ISFEs}
}
\references{
Bekhta, P., Ortynska, G. and Sedliacik, J. (2014). Properties of modified phenol-formaldehyde adhesive for plywood panels manufactured from high moisture content veneer. Drvna Industrija 65(4), 293-301.
}
\author{
Han Lin Shang
}
\keyword{methods}
