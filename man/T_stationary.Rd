\name{T_stationary}
\alias{T_stationary}
\title{
Testing stationarity of functional time series
}
\description{
A hypothesis test for stationarity of functional time series.
}
\usage{
T_stationary(sample, L = 49, J = 500, MC_rep = 1000, cumulative_var = .90, 
		    Ker1 = FALSE, Ker2 = TRUE, h = ncol(sample)^.5, pivotal = FALSE)
}
\arguments{
  \item{sample}{A discretised functional curves of dimension (p by n), where p represents the dimensionality and n represents sample size.}
  \item{L}{Number of Fourier basis functions.}
  \item{J}{Used to approximate the distribution of the squared integrals of Brownian bridges that appear in the limit distribution.}
  \item{MC_rep}{Number of replications.}
  \item{cumulative_var}{Amount of variance explained.}
  \item{Ker1}{Flat top kernel in (4.1) of Horvath et al. (2014).}
  \item{Ker2}{Flat top kernel in (7) of Politis (2003).}
  \item{h}{Bandwidth parameter.}
  \item{pivotal}{If \code{pivotal = TRUE}, a pivotal statistic is used.}
}
\details{
As in traditional (scalar and vector) time series analysis, the underlying assumption for inference in such models is stationarity. Stationarity is also required for functional dynamic regression models, for bootstrap and resampling methods for functional time series and for the functional analysis of volatility.
}
\value{
  \item{p_values}{When \code{p-value} is less than any level of significance, we reject the null hypothesis and conclude that the tested functional time series is not stationary.}
} 
\references{
L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional time series", \emph{Journal of Econometrics}, \bold{179}(1), 66-82.

D. N. Politis (2003) "Adaptive bandwidth choice", \emph{Journal of Nonparametric Statistics}, \bold{15}(4-5), 517-533.
}
\author{
Greg. Rice and Piotr. Kokoszka
}
\seealso{
\code{\link{farforecast}}
}
\examples{
N = 50
sampling_points = 10
sample = matrix(0,sampling_points+1,N)
for(ss in c(1:N))
{
     sample[,ss] = BBridge(x = 0, y = 0, t0 = 0, T = 1, N = sampling_points)
}
result = T_stationary(sample = sample, MC_rep=10)
}
\keyword{models}
