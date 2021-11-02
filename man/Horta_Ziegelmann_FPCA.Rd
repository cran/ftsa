\name{Horta_Ziegelmann_FPCA}
\alias{Horta_Ziegelmann_FPCA}
\title{
Dynamic functional principal component analysis for density forecasting
}
\description{
Implementation of a dynamic functional principal component analysis to forecast densities.  
}
\usage{
Horta_Ziegelmann_FPCA(data, gridpoints, h_scale = 1, p = 5, m = 5001, 
	kernel = c("gaussian", "epanechnikov"), band_choice = c("Silverman", "DPI"), 
	VAR_type = "both", lag_maximum = 6, no_boot = 1000, alpha_val = 0.1, 
	ncomp_select = "TRUE", D_val = 10)
}
\arguments{
  \item{data}{Densities or raw data matrix of dimension N by p, where N denotes sample size and p denotes dimensionality}
  \item{gridpoints}{Grid points}
  \item{h_scale}{Scaling parameter in the kernel density estimator}
  \item{p}{Number of backward parameters }
  \item{m}{Number of grid points}
  \item{kernel}{Type of kernel function}
  \item{band_choice}{Selection of optimal bandwidth}
  \item{VAR_type}{Type of vector autoregressive process}
  \item{lag_maximum}{A tuning parameter in the \code{super_fun} function}
  \item{no_boot}{A tuning parameter in the \code{super_fun} function}
  \item{alpha_val}{A tuning parameter in the \code{super_fun} function}
  \item{ncomp_select}{A tuning parameter in the \code{super_fun} function}
  \item{D_val}{A tuning parameter in the \code{super_fun} function}
}
\details{
1) Compute a kernel covariance function
2) Via eigen-decomposition, a density can be decomposed into a set of functional principal components and their associated scores
3) Fit a vector autoregressive model to the scores with the order selected by Akaike information criterion
4) By multiplying the estimated functional principal components with the forecast scores, obtain forecast densities
5) Since forecast densities may neither be non-negative nor sum to one, normalize the forecast densities accordingly
}
\value{
\item{Yhat.fix_den}{Forecast density}
\item{u}{Grid points}
\item{du}{Distance between two successive grid points}
\item{Ybar_est}{Mean of density functions}
\item{psihat_est}{Estimated functional principal components}
\item{etahat_est}{Estimated principal component scores}
\item{etahat_pred_val}{Forecast principal component scores}
\item{selected_d0}{Selected number of components}
\item{selected_d0_pvalues}{p-values associated with the selected functional principal components}
\item{thetahat_val}{Estimated eigenvalues}
}
\references{
Horta, E. and Ziegelmann, F. (2018) `Dynamics of financial returns densities: A functional approach applied to the Bovespa intraday index', \emph{International Journal of Forecasting}, \bold{34}, 75-88.

Bathia, N., Yao, Q. and Ziegelmann, F. (2010) `Identifying the finite dimensionality of curve time series', \emph{The Annals of Statistics}, \bold{38}, 3353-3386.
}
\author{
Han Lin Shang
}
\seealso{
\code{\link[ftsa]{CoDa_FPCA}}, \code{\link[ftsa]{LQDT_FPCA}}, \code{\link[ftsa]{skew_t_fun}}
}
\examples{
\dontrun{
Horta_Ziegelmann_FPCA(data = DJI_return, kernel = "epanechnikov", 
				band_choice = "DPI", ncomp_select = "FALSE")
}
}
\keyword{methods}
