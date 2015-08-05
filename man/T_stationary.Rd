\name{T_stationary}
\alias{T_stationary}
\title{
Testing stationarity of functional time series
}
\description{
A hypothesis test for stationarity of functional time series.
}
\usage{
T_stationary(sample, L = 49, J = 500, reps = 1000, seedlen = 50, 
		quan = c(0.9, 0.95), Ker1 = FALSE, Ker2 = TRUE)
}
\arguments{
  \item{sample}{A discretised functional curves of dimension (p by n), where p represents the dimensionality and n represents sample size}
  \item{L}{Number of Fourier basis functions}
  \item{J}{is used to approximate the distribution of the squared integrals of Brownian bridges that appear in the limit distribution.}
  \item{reps}{Number of replications}
  \item{seedlen}{It is used to define the "burn-in" length for the dependent processes that are used}
  \item{quan}{Selected quantiles}
  \item{Ker1}{Flat top kernel in (4.1) of Horvath et al. (2014)}
  \item{Ker2}{Flat top kernel in (7) of Politis (2003)}
}
\details{
As in traditional (scalar and vector) time series analysis, the underlying assumption for inference in such models is stationarity. Stationarity is also required for functional dynamic regression models, for bootstrap and resampling methods for functional time series and for the functional analysis of volatility.
}
\value{
  \item{reject9T}{Hypothesis testing results at 90\%. In the first column, the amount of the proportion of variation used for determining the number of components is 85\%. In the second column, the amount of the proportion of variation used for determining the number of components is 95\%. If a value is 0, this indicates non-stationarity of a functional time series for one repetition.}
  \item{reject95T}{Hypothesis testing results at 95\%. In the first column, the amount of the proportion of variation used for determining the number of components is 85\%. In the second column, the amount of the proportion of variation used for determining the number of components is 95\%. If a value is 0, this indicates non-stationarity of a functional time series for one repetition.}
} 
\references{
L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional time series", \emph{Journal of Econometrics}, \bold{179}, 66-82.

D. N. Politis (2003) "Adaptive bandwidth choice", \emph{Journal of Nonparametric Statistics}, \bold{25}, 517-533.
}
\author{
Greg. Rice, Piotr. Kokoszka and Han Lin Shang
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
dum = T_stationary(sample = sample, reps=10)
}
\keyword{models}
