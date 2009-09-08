\name{cm.spline}
\alias{cm.spline}
\title{Monotonic interpolating splines}
\description{
Performs cubic spline monotonic interpolation of given data points,
returning either a list of points obtained by the interpolation or a function performing the interpolation. The
splines are constrained to be monotonically increasing (i.e., the slope is never negative).}

\usage{
cm.spline(x, y = NULL, n = 3 * length(x), method = "fmm", xmin = min(x),
 xmax = max(x), gulim = 0)
}
\arguments{
\item{x,y}{Vectors giving the coordinates of the points to be interpolated. Alternatively a single plotting structure can be specified: see \code{\link[grDevices]{xy.coords}}.}
\item{method}{Specifies the type of spline to be used. Possible values are \dQuote{fmm}, \dQuote{natural} and \dQuote{periodic}.}
\item{n}{Interpolation takes place at \eqn{n} equally spaced points spanning the interval [\code{xmin}, \code{xmax}].}
\item{xmin}{Left-hand endpoint of the interpolation interval.}
\item{xmax}{Right-hand endpoint of the interpolation interval.}
\item{gulim}{If \code{gulim!=0} then it is taken as the upper limit on the gradient, and the interpolant is limited gradient rather than monotonic.}
}
\details{Results are identical to \code{\link[stats]{splinefun}} except that Hyman Filtering is used to produce co-monotonic
interpolation.
}

\value{
\item{spline}{Returns a list containing components \code{x} and \code{y} which give the coordinates where interpolation took place, and the interpolated values.}
\item{splinefun}{Returns a function which will perform cubic spline interpolation of the given data points. This is often more useful than \code{spline}.}
}

\references{
G. E. Forsythe and M. A. Malcolm and C. B. Moler (1977) \emph{Computer Methods for Mathematical Computations}, Englewood Cliffs, New Jersey: Prentice-Hall.

J. M. Hyman (1983) "Accurate monotonicity preserving cubic interpolation", \emph{SIAM Journal on Scientific Computing}, \bold{4}(4), 645-654.

R. L. Dougherty and A. Edelman and J. M. Hyman (1989) "Nonnegativity, monotonicity, or convexity-preserving cubic and quintic Hermit interpolation", \emph{Mathematics of Computation}, \bold{52}(186), 471-494.
}
\author{Simon N Wood, Rob J Hyndman}
\examples{
x <- seq(0, 4, l = 20)
y <- sort(rnorm(20))
plot(x,y)
lines(spline(x = x, y = y, n = 201), col = 2) # Not necessarily monotonic
lines(cm.spline(x = x, y = y, n = 201), col = 3) # Monotonic
}
\keyword{smooth}

