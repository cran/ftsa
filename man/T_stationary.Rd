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
		    Ker1 = FALSE, Ker2 = TRUE, h = ncol(sample)^.5, pivotal = FALSE,
		    use_table = FALSE, significance)
}
\arguments{
  \item{sample}{A matrix of discretised curves of dimension (p by n), where p represents the dimensionality and n represents sample size.}
  \item{L}{Number of Fourier basis functions.}
  \item{J}{Truncation level used to approximate the distribution of the squared integrals of Brownian bridges that appear in the limit distribution.}
  \item{MC_rep}{Number of replications.}
  \item{cumulative_var}{Amount of variance explained.}
  \item{Ker1}{Flat top kernel in (4.1) of Horvath et al. (2014).}
  \item{Ker2}{Flat top kernel in (7) of Politis (2003).}
  \item{h}{Kernel bandwidth.}
  \item{pivotal}{If \code{pivotal = TRUE}, a pivotal statistic is used.}
  \item{use_table}{If \code{use_table = TRUE}, use the critical values that are available in the book titled Inference for Functional Data (Table 6.1, page 88).}
  \item{significance}{Level of significance. Possibilities are \dQuote{10\%}, \dQuote{5\%}, \dQuote{1\%}.}
}
\details{
As in traditional (scalar and vector) time series analysis, many inferential procedures for functional time series assume stationarity. Stationarity is required for functional dynamic regression models, for bootstrap and resampling methods for functional time series and for the functional analysis of volatility.
}
\value{
  \item{p-value}{When \code{p-value} is less than any level of significance, we reject the null hypothesis and conclude that the tested functional time series is not stationary.}
}
\references{
L. Horvath and Kokoszka, P. (2012) Inference for Functional Data with Applications, Springer, New York.

L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional time series", \emph{Journal of Econometrics}, \bold{179}(1), 66-82.

D. N. Politis (2003) "Adaptive bandwidth choice", \emph{Journal of Nonparametric Statistics}, \bold{15}(4-5), 517-533.

A. Aue, G. Rice, O. S\"onmez (2018) "Detecting and dating structural breaks in functional data without dimension reduction", \emph{Journal of the Royal Statistical Society: Series B}, \bold{80}(3), 509-529.
}
\author{
Greg. Rice and Han Lin Shang
}
\seealso{
\code{\link{farforecast}}
}
\examples{
result = T_stationary(sample = pm_10_GR_sqrt$y)
result_pivotal = T_stationary(sample = pm_10_GR_sqrt$y, J = 100, MC_rep = 5000,
				h = 20, pivotal = TRUE)
}
\keyword{models}
