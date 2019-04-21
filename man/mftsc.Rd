\name{mftsc}
\alias{mftsc}
\title{
Clustering multiple functional time series
}
\description{
Clustering multiple functional time series. This function uses a functional panel data model to cluster different functional time series into subgroups.
}
\usage{
mftsc(X, alpha)
}
\arguments{
  \item{X}{A list of functional time series to be clustered for each object, it is a \code{p} by \code{q} matrix, where \code{p} denotes sample size and \code{q} denotes number of grid points within a function support.}
  \item{alpha}{An input value for adjusted rand index to measure similarity of the memberships with last iteration, can be any value bigger than 0.9.}
}
\value{
\item{iteration}{Number of iterations until convergence}
\item{membership}{A list of all membership matrices at each iteration}
\item{member.final}{Final membership}
}
\author{
Chen Tang and Han Lin Shang
}
\seealso{
\code{\link{dmfpca}}
}
\examples{
## The following takes about 20 minutes to run ##
\dontrun{
cluster_result <- mftsc(X = sim_ex_cluster, alpha = 0.99)
}
}
\keyword{methods}
