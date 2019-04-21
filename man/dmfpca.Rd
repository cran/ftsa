\name{dmfpca}
\alias{dmfpca}
\title{
Dynamic multilevel functional principal component analysis
}
\description{
Functional principal component analysis is used to decompose multiple functional time series. This function uses a functional panel data model to reduce dimensions for multiple functional time series.
}
\usage{
dmfpca(y, M = NULL, J = NULL, N = NULL, tstart = 0, tlength = 1)
}
\arguments{
  \item{y}{A data matrix containing functional responses. Each row contains measurements from a function at a set of grid points, and each column contains measurements of all functions at a particular grid point}
  \item{M}{Number of \code{fts} obejcts}
  \item{J}{Number of functions in each object}
  \item{N}{Number of grid points per function}
  \item{tstart}{Start point of the grid points}
  \item{tlength}{Length of the interval that the functions are evaluated at}
}
\value{
\item{K1}{Number of components for the common time-trend}
\item{K2}{Number of components for the residual component}
\item{lambda1}{A vector containing all common time-trend eigenvalues in non-increasing order}
\item{lambda2}{A vector containing all residual component eigenvalues in non-increasing order}
\item{phi1}{A matrix containing all common time-trend eigenfunctions. Each row contains an eigenfunction evaluated at the same set of grid points as the input data. The eigenfunctions are in the same order as the corresponding eigenvalues}
\item{phi2}{A matrix containing all residual component eigenfunctions. Each row contains an eigenfunction               evaluated at the same set of grid points as the input data.  The eigenfunctions are in the same order as the corresponding eigenvalues.}
\item{scores1}{A matrix containing estimated common time-trend principal component scores. Each row corresponding to the common time-trend scores for a particular subject in a cluster. The number of rows is the same as that of the input matrix y. Each column contains the scores for a common time-trend component for all subjects.}
\item{scores2}{A matrix containing estimated residual component principal component scores. Each row corresponding to the level 2 scores for a particular subject in a cluster.  The number of rows is the same as that of the input matrix y. Each column contains the scores for a residual component for all subjects.}
\item{mu}{A vector containing the overall mean function.}
\item{eta}{A matrix containing the deviation from overall mean function to country specific mean function. The number of rows is the number of countries.}
}
\references{
Rice, G. and Shang, H. L. (2017) "A plug-in bandwidth selection procedure for long-run covariance estimation with stationary functional time series", \emph{Journal of Time Series Analysis}, \bold{38}, 591-609.

Shang, H. L. (2016) "Mortality and life expectancy forecasting for a group of populations in developed countries: A multilevel functional data method", \emph{The Annals of Applied Statistics}, \bold{10}, 1639-1672.

Di, C.-Z., Crainiceanu, C. M., Caffo, B. S. and Punjabi, N. M. (2009) "Multilevel functional principal component analysis", \emph{The Annals of Applied Statistics}, \bold{3}, 458-488.
}
\author{
Chen Tang and Han Lin Shang
}
\seealso{
\code{\link{mftsc}}
}
\examples{
## The following takes about 10 seconds to run ##
\dontrun{
y <- do.call(rbind, sim_ex_cluster) 
MFPCA.sim <- dmfpca(y, M = length(sim_ex_cluster), J = nrow(sim_ex_cluster[[1]]), 
				    N = ncol(sim_ex_cluster[[1]]), tlength = 1)
}
}
\keyword{methods}
