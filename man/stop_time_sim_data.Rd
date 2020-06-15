\name{stop_time_sim_data}
\alias{stop_time_sim_data}
\title{
Simulated functional time series from a functional autoregression of order one
}
\description{
For detecting the optimal stopping time, we simulate a curve time series that follows a functional autoregression of order 1, with a breakpoint in the middle point of the entire sample.
}
\usage{
stop_time_sim_data(sample_size, omega, seed_number)
}
\arguments{
  \item{sample_size}{Number of curves}
  \item{omega}{Noise level}
  \item{seed_number}{Random seed number}
}
\value{
An object of class \code{fts}
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{stop_time_detect}}
}
\examples{
stop_time_sim_data(sample_size = 401, omega = 0.1, seed_number = 123)
}
\keyword{datasets}
