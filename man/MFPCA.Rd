\name{MFPCA}
\alias{MFPCA}
\title{
Multilevel functional principal component analysis for clustering
}
\description{
A multilevel functional principal component analysis for performing clustering analysis
}
\usage{
MFPCA(y, M = NULL, J = NULL, N = NULL)
}
\arguments{
  \item{y}{A data matrix containing functional responses. Each row contains measurements from a function at a set of grid points, and each column contains measurements of all functions at a particular grid point}
  \item{M}{Number of countries}
  \item{J}{Number of functional responses in each country}
  \item{N}{Number of grid points per function}
}
\value{
\item{K1}{Number of components at level 1}
\item{K2}{Number of components at level 2}
\item{K3}{Number of components at level 3}
\item{lambda1}{A vector containing all level 1 eigenvalues in non-increasing order}
\item{lambda2}{A vector containing all level 2 eigenvalues in non-increasing order}
\item{lambda3}{A vector containing all level 3 eigenvalues in non-increasing order}
\item{phi1}{A matrix containing all level 1 eigenfunctions. Each row contains an eigenfunction evaluated at the same set of grid points as the input data. The eigenfunctions are in the same order as the corresponding eigenvalues}
\item{phi2}{A matrix containing all level 2 eigenfunctions. Each row contains an eigenfunction evaluated at the same set of grid points as the input data. The eigenfunctions are in the same order as the corresponding eigenvalues}
\item{phi3}{A matrix containing all level 3 eigenfunctions. Each row contains an eigenfunction evaluated at the same set of grid points as the input data. The eigenfunctions are in the same order as the corresponding eigenvalues}
\item{scores1}{A matrix containing estimated level 1 principal component scores. Each row corresponds to the level 1 scores for a particular subject in a cluster. The number of rows is the same as that of the input matrix \code{y}. Each column contains the scores for a level 1 component for all subjects}
\item{scores2}{A matrix containing estimated level 2 principal component scores. Each row corresponds to the level 2 scores for a particular subject in a cluster. The number of rows is the same as that of the input matrix \code{y}. Each column contains the scores for a level 2 component for all subjects.}
\item{scores3}{A matrix containing estimated level 3 principal component scores. Each row corresponds to the level 3 scores for a particular subject in a cluster. The number of rows is the same as that of the input matrix \code{y}. Each column contains the scores for a level 3 component for all subjects.}
\item{mu}{A vector containing the overall mean function}
\item{eta}{A matrix containing the deviation from overall mean function to country-specific mean function. The number of rows is the number of countries}
\item{Rj}{Common trend}
\item{Uij}{Country-specific mean function}
}
\author{
Chen Tang, Yanrong Yang and Han Lin Shang
}
\seealso{
\code{\link{mftsc}}
}
\keyword{methods}
