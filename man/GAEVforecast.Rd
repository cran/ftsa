\name{GAEVforecast}
\alias{GAEVforecast}
\title{
Fit a generalized additive extreme value model to the functional data with given basis numbers
}
\description{
One-step-ahead forecast for any given quantile(s) of functional time sereies of extreme values using a generalized additive extreme value (GAEV) model.
}
\usage{
GAEVforecast(data, q, d.loc.max = 10, d.logscale.max = 10)
}
\arguments{
  \item{data}{a n by p data matrix, where n denotes the number of functional objects and p denotes the number of realizations on each functional object}
  \item{q}{a required scalar or vector of GEV quantiles that are of forecasting interest}
  \item{d.loc.max}{the maximum number of basis functions considered for the location parameter}
  \item{d.logscale.max}{the maximum number of basis functions considered for the (log-)scale parameter}
}
\details{
For the functional time seres \eqn{\{X_t(u),t=1,...,T,u\in \mathcal{I}\}}, the GAEV model is given as
\deqn{
X_{t}(u) ~ GEV[\mu_{t}(u),\sigma_t(u),\xi_t], 
}
where
\deqn{
\mu_t(u) = \beta^{(\mu)}_{t,0} + \sum_{i=1}^{d_1}\beta^{(\mu)}_{t,i}b^{(\mu)}_{i}(u),
}
\deqn{
\ln(\sigma_t(u)) = \beta^{(\sigma)}_{t,0} + \sum_{i=1}^{d_2}\beta^{(\sigma)}_{t,i}b^{(\sigma)}_{i}(u), \xi_t \in [0,\infty),
}
where \eqn{d_{j},j=1,2} are positive integers of basis numbers, \eqn{\{b^{(\mu)}_{i}(u),i=1,\dots,d_{1}\}} and \eqn{\{b^{(\sigma)}_{i}(u),i=1,\dots,d_{2}\}} are the cubic regression spline basis functions.

The optimal number of basis functions \eqn{(d_1,d_2)} are chosen by minimizing the Kullback-Leibler divergence on the test set using a leave-one-out cross-validation technique.

The one-step-ahead forecast of the joint coefficients \eqn{(\widehat{\beta^{(\mu)}}_{T+1,i},\widehat{\beta^{(\sigma)}}_{T+1,j},\widehat{\xi}_{T+1},i=0,...,d_1,j=0,...,d_2)} are produced using a vector autoregressive model, whose order is selected via the corrected Akaike information criterion. Then the one-step-ahead forecast of the GEV parameter \eqn{(\widehat{\mu}_{T+1}(u),\widehat{\sigma}_{T+1}(u),\widehat{\xi}_{T+1})} can be computed accordingly.

The one-step-ahead forecast for the \eqn{\tau}-th quantile of the extreme values \eqn{\widehat{X}_{T+1}(u)} is computed by
\deqn{
    Q_{\tau}(u|\widehat{\mu}_{T+1},\widehat{\sigma}_{T+1},\widehat{\xi}_{T+1})}
   = \deqn{
   \widehat{\mu}_{T+1}(u) + \frac{\widehat{\sigma}_{T+1}(u) \big[(-\ln(\tau))^{-\widehat{\xi}_{T+1}}-1\big]}{\widehat{\xi}_{T+1}},  \xi > 0, \tau\in [0,1);\ \xi < 0, \tau\in (0,1], \\
   \widehat{\mu}_{T+1}(u) - \widehat{\sigma}_{T+1}(u) \cdot \ln[-\ln\big(\tau)], \xi=0, \tau \in (0,1).
}
}
\value{
\item{kdf.location}{the optimal number of basis functions considered for the location parameter}
\item{kdf.logscale}{the optimal number of basis functions considered for the (log-)scale parameter}
\item{basis.location}{the basis functions for the location parameter}
\item{basis.logscale}{the basis functions for the (log-)scale parameter}
\item{para.location.pred}{the predicted location function}
\item{para.scale.pred}{the predicted scale function}
\item{para.shape.pred}{the predicted shape parameter}
\item{density.pred}{the prediced density function(s) for the given quantile(s)}
}
\references{
Shang, H. L. and Xu, R. (2021) `Functional time series forecasting of extreme values', \emph{Communications in Statistics Case Studies Data Analysis and Applications}, \bold{in press}.
}
\author{
Ruofan Xu and Han Lin Shang
}
\examples{
\dontrun{
library(evd)
data = matrix(rgev(1000),ncol=50) 
GAEVforecast(data = data, q = c(0.02,0.7), d.loc.max = 5, d.logscale.max = 5)
}
}
\keyword{methods}
