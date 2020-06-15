\name{MAF_multivariate}
\alias{MAF_multivariate}
\title{
Maximum autocorrelation factors
}
\description{
Dimension reduction via maximum autocorrelation factors
}
\usage{
MAF_multivariate(data, threshold)
}
\arguments{
  \item{data}{A p by n data matrix, where p denotes the number of variables and n denotes the sample size}
  \item{threshold}{A threshold level for retaining the optimal number of factors}
}
\value{
\item{MAF}{Maximum autocorrelation factor scores}
\item{MAF_loading}{Maximum autocorrelation factors}
\item{Z}{Standardized original data}
\item{recon}{Reconstruction via maximum autocorrelation factors}
\item{recon_err}{Reconstruction errors between the standardized original data and reconstruction via maximum autocorrelation factors}
\item{ncomp_threshold}{Number of maximum autocorrelation factors selected by explaining autocorrelation at and above a given level of threshold}
\item{ncomp_eigen_ratio}{Number of maximum autocorrelation factors selected by eigenvalue ratio tests}
}
\references{
M. A. Haugen, B. Rajaratnam and P. Switzer (2015). Extracting common time trends from concurrent time series: Maximum autocorrelation factors with applications, arXiv paper https://arxiv.org/abs/1502.01073.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link{ftsm}}
}
\examples{
MAF_multivariate(data = pm_10_GR_sqrt$y, threshold = 0.85)
}
\keyword{methods}
