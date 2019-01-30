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
A. Cuevas and M. Febrero and R. Fraiman (2006), "On the use of the bootstrap for estimating functions with functional data", \emph{Computational Statistics and Data Analysis}, \bold{51}(2), 1063-1074.

A. Cuevas and M. Febrero and R. Fraiman (2007), "Robust estimation and classification for functional data via projection-based depth notions", \emph{Computational Statistics}, \bold{22}(3), 481-496.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2007) "A functional analysis of NOx levels: location and scale estimation and outlier detection", \emph{Computational Statistics}, \bold{22}(3), 411-427.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2008) "Outlier detection in functional data by depth measures, with application to identify abnormal NOx levels", \emph{Environmetrics}, \bold{19}(4), 331-345.

M. Febrero and P. Galeano and W. Gonzalez-Manteiga (2010) "Measures of influence for the functional linear model with scalar response", \emph{Journal of Multivariate Analysis}, \bold{101}(2), 327-339. 

J. A. Cuesta-Albertos and A. Nieto-Reyes (2010) "Functional classification and the random Tukey depth. Practical issues", Combining Soft Computing and Statistical Methods in Data Analysis, \emph{Advances in Intelligent and Soft Computing}, \bold{77}, 123-130.

D. Gervini (2012) "Outlier detection and trimmed estimation in general functional spaces", \emph{Statistica Sinica}, \bold{22}(4), 1639-1660.

H. L. Shang (2015) "Re-sampling techniques for estimating the distribution of descriptive statistics of functional data", \emph{Communication in Statistics--Simulation and Computation}, \bold{44}(3), 614-635.

H. L. Shang (2018) Bootstrap methods for stationary functional time series, \emph{Statistics and Computing}, \bold{28}(1), 1-10.
}
\author{Han Lin Shang}
\seealso{
\code{\link[ftsa]{pcscorebootstrapdata}}
}
\examples{
# Bootstrapping the distribution of a summary statistics of functional data.
fbootstrap(data = ElNino_ERSST_region_1and2)
}
\keyword{multivariate}

