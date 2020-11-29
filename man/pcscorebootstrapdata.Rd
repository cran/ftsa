\name{pcscorebootstrapdata}
\alias{pcscorebootstrapdata}
\title{
Bootstrap independent and identically distributed functional data or functional time series
}
\description{
Computes bootstrap or smoothed bootstrap samples based on either independent and identically distributed functional data or functional time series.
}
\usage{
pcscorebootstrapdata(dat, bootrep, statistic, bootmethod = c("st", "sm", 
	"mvn", "stiefel", "meboot"), smo)
}
\arguments{
  \item{dat}{An object of class \code{matrix}.}
  \item{bootrep}{Number of bootstrap samples.}
  \item{statistic}{Summary statistics.}
  \item{bootmethod}{Bootstrap method. When \code{bootmethod = "st"}, the sampling with replacement is implemented. To avoid the repeated bootstrap samples, 
  the smoothed boostrap method can be implemented by adding multivariate Gaussian random noise. When \code{bootmethod = "mvn"}, the bootstrapped principal 
  component scores are drawn from a multivariate Gaussian distribution with the mean and covariance matrices of the original principal component scores. 
  When \code{bootmethod = "stiefel"}, the bootstrapped principal component scores are drawn from a Stiefel manifold with the mean and covariance matrices of
  the original principal component scores. When \code{bootmethod = "meboot"}, the bootstrapped principal component scores are drawn from a maximum entropy algorithm of Vinod (2004).}
  \item{smo}{Smoothing parameter.}
}
\details{
We will presume that each curve is observed on a grid of \eqn{T} points with \eqn{0\leq t_1<t_2\dots<t_T\leq \tau}. 
Thus, the raw data set \eqn{(X_1,X_2,\dots,X_n)} of \eqn{n} observations will consist of an \eqn{n} by \eqn{T} data matrix.
By applying the singular value decomposition, \eqn{X_1,X_2,\dots,X_n} can be decomposed into \eqn{X = ULR^{\top}}, 
where the crossproduct of \eqn{U} and \eqn{R} is identity matrix.  

Holding the mean and \eqn{L} and \eqn{R} fixed at their realized values, there are four re-sampling methods that differ mainly by the ways of re-sampling U.

(a) Obtain the re-sampled singular column matrix by randomly sampling with replacement from the original principal component scores.

(b) To avoid the appearance of repeated values in bootstrapped principal component scores, we adapt a smooth bootstrap procedure by adding a white noise component to the bootstrap.

(c) Because principal component scores follow a standard multivariate normal distribution asymptotically, we can randomly draw principal component scores from a multivariate normal distribution with mean vector and covariance matrix of original principal component scores.

(d) Because the crossproduct of U is identitiy matrix, U is considered as a point on the Stiefel manifold, that is the space of \eqn{n} orthogonal vectors, thus we can randomly draw principal component scores from the Stiefel manifold.
}
\value{
\item{bootdata}{Bootstrap samples. If the original data matrix is \eqn{p} by \eqn{n}, then the bootstrapped data are \eqn{p} by \eqn{n} by \eqn{bootrep}.}
\item{meanfunction}{Bootstrap summary statistics. If the original data matrix is \eqn{p} by \eqn{n}, then the bootstrapped summary statistics is \eqn{p} by \eqn{bootrep}.}
}
\references{
H. D. Vinod (2004), "Ranking mutual funds using unconventional utility theory and stochastic dominance", \emph{Journal of Empirical Finance}, \bold{11}(3), 353-377.

A. Cuevas, M. Febrero, R. Fraiman (2006), "On the use of the bootstrap for estimating functions with functional data", \emph{Computational Statistics and Data Analysis}, \bold{51}(2), 1063-1074.

D. S. Poskitt and A. Sengarapillai (2013), "Description length and dimensionality reduction in functional data analysis", \emph{Computational Statistics and Data Analysis}, \bold{58}, 98-113.

H. L. Shang (2015), "Re-sampling techniques for estimating the distribution of descriptive statistics of functional data", \emph{Communications in Statistics--Simulation and Computation}, \bold{44}(3), 614-635.

H. L. Shang (2018), "Bootstrap methods for stationary functional time series", \emph{Statistics and Computing}, \bold{28}(1), 1-10.
}
\author{
Han Lin Shang
}

\seealso{
\code{\link[ftsa]{fbootstrap}}
}
\examples{
# Bootstrapping the distribution of a summary statistics of functional data.	
boot1 = pcscorebootstrapdata(dat = ElNino_ERSST_region_1and2$y, bootrep = 200, 
	statistic = "mean", bootmethod = "st")
boot2 = pcscorebootstrapdata(dat = ElNino_ERSST_region_1and2$y, bootrep = 200, 
	statistic = "mean", bootmethod = "sm", smo = 0.05)
boot3 = pcscorebootstrapdata(dat = ElNino_ERSST_region_1and2$y, bootrep = 200, 
	statistic = "mean", bootmethod = "mvn")
boot4 = pcscorebootstrapdata(dat = ElNino_ERSST_region_1and2$y, bootrep = 200, 
	statistic = "mean", bootmethod = "stiefel")
boot5 = pcscorebootstrapdata(dat = ElNino_ERSST_region_1and2$y, bootrep = 200, 
	statistic = "mean", bootmethod = "meboot")
}
\keyword{models}
