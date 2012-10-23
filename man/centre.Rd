\name{centre}
\alias{centre}
\title{
Mean function, variance function, median function, trim mean function of functional data
}
\description{
Mean function, variance function, median function, trim mean function of functional data
}
\usage{
centre(x, type)
}
\arguments{
  \item{x}{An object of class \code{matrix}.}
  \item{type}{Mean, variance, median or trim mean?}
}
\value{
Return mean function, variance function, median function or trim mean function.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{pcscorebootstrapdata}}
}
\examples{
# mean function is often removed in the functional principal component analysis.
# trimmed mean function is sometimes employed for robustness in the presence of outliers.
# In calculating trimmed mean function, several functional depth measures were employed.	
centre(ElNino$y, "mean")
}
\keyword{models}

