\name{LQDT_FPCA}
\alias{LQDT_FPCA}
\title{
Log quantile density transform
}
\description{
Probability density function, cumulative distribution function and quantile density function are three characterizations of a distribution. Of these three, quantile density function is the least constrained. The only constrain is nonnegative. By taking a log transformation, there is no constrain. 
}
\usage{
LQDT_FPCA(data, gridpoints, h_scale = 1, M = 3001, m = 5001, lag_maximum = 4, 
		no_boot = 1000, alpha_val = 0.1, p = 5, 
		band_choice = c("Silverman", "DPI"), 
		kernel = c("gaussian", "epanechnikov"), 
		forecasting_method = c("uni", "multi"), 
		varprop = 0.85, fmethod, VAR_type)
}
\arguments{
  \item{data}{Densities or raw data matrix of dimension N by p, where N denotes sample size and p denotes dimensionality}
  \item{gridpoints}{Grid points}
  \item{h_scale}{Scaling parameter in the kernel density estimator}
  \item{M}{Number of grid points between 0 and 1}
  \item{m}{Number of grid points within the data range}
  \item{lag_maximum}{A tuning parameter in the \code{super_fun} function}
  \item{no_boot}{A tuning parameter in the \code{super_fun} function}
  \item{alpha_val}{A tuning parameter in the \code{super_fun} function}
  \item{p}{Number of backward parameters}
  \item{band_choice}{Selection of optimal bandwidth}
  \item{kernel}{Type of kernel function}
  \item{forecasting_method}{Univariate or multivariate time series forecasting method}
  \item{varprop}{Proportion of variance explained}
  \item{fmethod}{If \code{forecasting_method = "uni"}, specify a particular forecasting method}
  \item{VAR_type}{If \code{forecasting_method = "multi"}, specify a particular type of vector autoregressive model}
}
\details{
1) Transform the densities f into log quantile densities Y and c specifying the value of the cdf at 0 for the target density f.
2) Compute the predictions for future log quantile density and c value.
3) Transform the forecasts in Step 2) into the predicted density f.
}
\value{
\item{L2Diff}{L2 norm difference between reconstructed and actual densities}
\item{unifDiff}{Uniform Metric excluding missing boundary values (due to boundary cutoff)}
\item{density_reconstruct}{Reconstructed densities}
\item{density_original}{Actual densities}
\item{dens_fore}{Forecast densities}
\item{totalMass}{Assess loss of mass incurred by boundary cutoff}
\item{u}{m number of grid points} 
}
\references{
Petersen, A. and Muller, H.-G. (2016) `Functional data analysis for density functions by transformation to a Hilbert space', \emph{The Annals of Statistics}, \bold{44}, 183-218.

Jones, M. C. (1992) `Estimating densities, quantiles, quantile densities and density quantiles', \emph{Annals of the Institute of Statistical Mathematics}, \bold{44}, 721-727.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{CoDa_FPCA}}, \code{\link[ftsa]{Horta_Ziegelmann_FPCA}}, \code{\link[ftsa]{skew_t_fun}}
}
\examples{
\dontrun{
LQDT_FPCA(data = DJI_return, band_choice = "DPI", kernel = "epanechnikov", 
			forecasting_method = "uni", fmethod = "ets")
}		
}
\keyword{methods}
