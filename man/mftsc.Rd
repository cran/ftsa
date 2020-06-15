\name{mftsc}
\alias{mftsc}
\title{
 Multiple funtional time series clustering
}
\description{
Clustering the multiple functional time series. The function uses the functional panel data model to cluster different time series into subgroups
}
\usage{mftsc(X, alpha)}
\arguments{
\item{X}{A list of sets of smoothed functional time series to be clustered, for each object, it is a p x q matrix, where p is the sample size and q is the number of grid points of the function}

\item{alpha}{A value input for adjusted rand index to measure similarity of the memberships with last iteration, can be any value big than 0.9}
}
\details{As an initial step, conventional k-means clustering is performed on the dynamic FPC scores, then an iterative membership updating process is applied by fitting the MFPCA model.}
\value{
\item{iteration}{the number of iterations until convergence}
\item{memebership}{a list of all the membership matrices at each iteration}
\item{member.final}{the final membership}
}
\author{Chen Tang, Yanrong Yang and Han Lin Shang}
\seealso{
\code{\link{MFPCA}}
}
\examples{
\dontrun{
data(sim_ex_cluster)
cluster_result<-mftsc(X=sim_ex_cluster, alpha=0.99)
cluster_result$member.final
}
}
