\name{fbootstrap}
\alias{fbootstrap}
\title{Bootstrap independent and identically distributed functional data}
\description{
Computes bootstrap or smoothed bootstrap samples based on independent and identically distributed functional data.
}
\usage{
fbootstrap(data, estad = func.mean, alpha = 0.05, nb = 200, suav = 0,
 media.dist = FALSE, graph = FALSE, ...)
}
\arguments{
  \item{data}{An object of class \code{fds} or \code{fts}.}
  \item{estad}{Estimate function of interest. Default is to estimate the mean function. Other options are \code{func.mode} or \code{func.var}.}
  \item{alpha}{Significance level used in the smooth bootstrapping.}
  \item{nb}{Number of bootstrap samples.}
  \item{suav}{Smoothing parameter.}
  \item{media.dist}{Estimate mean function.}
  \item{graph}{Graphical output.}
  \item{\dots}{Other arguments.}
}
\value{
A list containing the following components is returned.
  \item{estimate}{Estimate function.}
  \item{max.dist}{Max distance of bootstrap samples.}
  \item{rep.dist}{Distances of bootstrap samples.}
  \item{resamples}{Bootstrap samples.}
  \item{center}{Functional mean.}
}
\references{
A. Cuevas and M. Febrero and R. Fraiman (2006), "On the use of the bootstrap for estimating functions with functional data", Computational Statistics and Data Analysis, \bold{51}(2), 1063-1074.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2007) "A functional analysis of NOx levels: location and scale estimation and outlier detection", \emph{Computational Statistics}, \bold{22}(3), 411-427.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2010) "Measures of influence for the functional linear model with scalar response", \emph{Journal of Multivariate Analysis}, \bold{101}(2), 327-339. 
}
\author{Han Lin Shang}
\seealso{
\code{\link[ftsa]{pcscorebootstrapdata}}
}
\examples{
fbootstrap(data = ElNino)
}
\keyword{multivariate}

