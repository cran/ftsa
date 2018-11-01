\name{MFDM}
\alias{MFDM}
\title{
Multilevel functional data method
}
\description{
Fit a multilevel functional principal component model. The function uses two-step functional principal component decompositions.
}
\usage{
MFDM(mort_female, mort_male, mort_ave, percent_1 = 0.95, percent_2 = 0.95, fh, 
	     level = 80, alpha = 0.2, MCMCiter = 100, fmethod = c("auto_arima", "ets"), 
		   BC = c(FALSE, TRUE), lambda)
}
\arguments{
  \item{mort_female}{Female mortality (p by n matrix), where p denotes the dimension and n denotes the sample size.}
  \item{mort_male}{Male mortality (p by n matrix).}
  \item{mort_ave}{Total mortality (p by n matrix).}
  \item{percent_1}{Cumulative percentage used for determining the number of common functional principal components.}
  \item{percent_2}{Cumulative percentage used for determining the number of sex-specific functional principal components.}
  \item{fh}{Forecast horizon.}
  \item{level}{Nominal coverage probability of a prediction interval.}
  \item{alpha}{1 - Nominal coverage probability.}
  \item{MCMCiter}{Number of MCMC iterations.}
  \item{fmethod}{Univariate time-series forecasting method.}
  \item{BC}{If Box-Cox transformation is performed.}
  \item{lambda}{If \code{BC = TRUE}, specify a Box-Cox transformation parameter.}
}
\details{
The basic idea of multilevel functional data method is to decompose functions from different sub-populations into an aggregated average, a sex-specific deviation from the aggregated average, a common trend, a sex-specific trend and measurement error. The common and sex-specific trends are modelled by projecting them onto the eigenvectors of covariance operators of the aggregated and sex-specific centred stochastic process, respectively.
}
\value{
\item{first_percent}{Percentage of total variation explained by the first common functional principal component.}
\item{female_percent}{Percentage of total variation explained by the first female functional principal component in the residual.}
\item{male_percent}{Percentage of total variation explained by the first male functional principal component in the residual.}
\item{mort_female_fore}{Forecast female mortality in the original scale.}
\item{mort_male_fore}{Forecast male mortality in the original scale.}
}
\references{
C. M. Crainiceanu and J. Goldsmith (2010) "Bayesian functional data analysis using WinBUGS", \emph{Journal of Statistical Software}, \bold{32}(11).

C-Z. Di and C. M. Crainiceanu and B. S. Caffo and N. M. Punjabi (2009) "Multilevel functional principal component analysis", \emph{The Annals of Applied Statistics}, \bold{3}(1), 458-488.

V. Zipunnikov and B. Caffo and D. M. Yousem and C. Davatzikos and B. S. Schwartz and C. Crainiceanu (2015) "Multilevel functional principal component analysis for high-dimensional data", \emph{Journal of Computational and Graphical Statistics}, \bold{20}, 852-873.

H. L. Shang, P. W. F. Smith, J. Bijak, A. Wisniowski (2016) "A multilevel functional data method for forecasting population, with an application to the United Kingdom", \emph{International Journal of Forecasting}, \bold{32}(3), 629-649.
}
\author{
Han Lin Shang
}
\note{
It can be quite time consuming, especially when MCMCiter is large.
}
\seealso{
\code{\link{ftsm}}, \code{\link{forecast.ftsm}}, \code{\link{fplsr}}, \code{\link{forecastfplsr}}
}
\keyword{models}
