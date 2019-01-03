\name{hdfpca}
\alias{hdfpca}
\title{
High-dimensional functional principal component analysis
}
\description{
Fit a high dimensional functional principal component analysis model to a multiple-population of functional time series data.
}
\usage{
hdfpca(y, order, q = sqrt(dim(y[[1]])[2]), r)
}
\arguments{
  \item{y}{A list, where each item is a population of functional time series. Each item is a data matrix of dimension p by n, where p is the number of discrete points in each function and n is the sample size}
  \item{order}{The number of principal component scores to retain in the first step dimension reduction}
  \item{q}{The tuning parameter used in the first step dimension reduction, by default it is equal to the square root of the sample size}
  \item{r}{The number of factors to retain in the second step dimension reduction}
}
\details{
In the first step, dynamic functional principal component analysis is performed on each population and then in the second step, factor models are fitted to the resulting principal component scores. The high-dimensional functional time series are thus reduced to low-dimensional factors.
}
\value{
\item{y}{The input data}
\item{p}{The number of discrete points in each function}
\item{fitted}{A list containing the fitted functions for each population}
\item{m}{The number of populations}
\item{model}{Model 1 includes the first step dynamic functional principal component analysis models, model 2 includes the second step high-dimensional principal component analysis models}
\item{order}{Input order}
\item{r}{Input r}
}
\references{
Y. Gao, H. L. Shang and Y. Yang (2018) High-dimensional functional time series forecasting: An application to age-specific mortality rates, \emph{Journal of Multivariate Analysis}, \bold{forthcoming}.
}
\author{
Y. Gao and H. L. Shang
}
\seealso{
\code{\link{forecast.hdfpca}}, \code{\link{hd_data}}
}
\examples{
hd_model = hdfpca(hd_data, order = 2, r = 2)
}
\keyword{models}