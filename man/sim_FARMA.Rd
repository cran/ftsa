\name{sim_FARMA}
\alias{sim_FARMA}
\title{
Simulating data from a functional autoregressive moving average process
}
\description{
Simulating data from a stationary functional autoregressive moving average process
}
\usage{
sim_FARMA(n = 500, d = 0.1, no_grid = 101, seed_number, process = c("FAR", "FARMA"))
}
\arguments{
  \item{n}{Number of curves}
  \item{d}{Memory parameter}
  \item{no_grid}{Number of grid points within a unit interval}
  \item{seed_number}{Pseudo-random seed number}
  \item{process}{FAR or FARMA process}
}
\value{
Simulated curves
}
\references{
D. Li, P. M. Robinson, H. L. Shang (2020) Long-range dependent curve time series. Journal of the American Statistical Association, 115(530), 957-971.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{H_estimate_RS}}
}
\examples{
sim_data_set = sim_FARMA(seed_number = 123, process = "FAR")
}
\keyword{time series}
