\name{ER_GR}
\alias{ER_GR}
\title{
Selection of the number of principal components
}
\description{
Eigenvalue ratio and growth ratio
}
\usage{
ER_GR(data)
}
\arguments{
  \item{data}{An n by p matrix, where n denotes sample size and p denotes the number of discretized data points in a curve}
}
\value{
\item{k_ER}{The number of components selected by the eigenvalue ratio}
\item{k_GR}{The number of components selected by the growth ratio}
}
\references{
Lam, C. and Yao, Q. (2012). Factor modelling for high-dimensional time series: Inference for the number of factors. The Annals of Statistics, 40, 694-726.

Ahn, S. and Horenstein, A. (2013). Eigenvalue ratio test for the number of factors. Econometrica, 81, 1203-1227.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{ftsm}}
}
\examples{
ER_GR(pm_10_GR$y)
}
\keyword{methods}
