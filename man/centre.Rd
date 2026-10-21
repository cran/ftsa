\name{centre}
\alias{centre}
\title{
Mean function, variance function, median function, trimmed mean function of functional data
}
\description{
Mean function, variance function, median function, trimmed mean function of functional data
}
\usage{
centre(x, type)
}
\arguments{
  \item{x}{An object of class \code{matrix}.}
  \item{type}{Mean, variance, median, or trimmed mean?}
}
\value{
Return mean function, variance function, median function, or trimmed mean function.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{pcscorebootstrapdata}}, \code{\link[ftsa]{mean.fts}}, \code{\link[ftsa]{median.fts}}, \code{\link[ftsa]{sd.fts}}, \code{\link[ftsa]{var.fts}}
}
\examples{
# The mean function is often removed in functional principal component analysis.
# The trimmed mean function is sometimes employed for robustness in the presence of outliers.
# In calculating the trimmed mean function, several functional depth measures were employed.	
centre(x = ElNino_ERSST_region_1and2$y, type = "mean")
centre(x = ElNino_ERSST_region_1and2$y, type = "var")
centre(x = ElNino_ERSST_region_1and2$y, type = "median")
centre(x = ElNino_ERSST_region_1and2$y, type = "trimmed")
}
\keyword{models}

