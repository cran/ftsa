\name{fplsr}
\alias{fplsr}
\title{Functional partial least squares regression}
\description{
Fits a functional partial least squares (PLSR) model using nonlinear partial least squares
(NIPALS) algorithm or simple partial least squares (SIMPLS) algorithm.
}
\usage{
fplsr(data, order = 6, type = c("simpls", "nipals"), unit.weights =
 TRUE, weight = FALSE, beta = 0.1, interval = FALSE, method =
  c("delta", "boota"), alpha = 0.05, B = 100, adjust = FALSE,
   backh = 10)
}
\arguments{
  \item{data}{An object of class \code{fts}.}
  \item{order}{Number of principal components to fit.}
  \item{type}{When \code{type = "nipals"}, uses the NIPALS algorithm; when \code{type = "simpls"}, uses the SIMPLS algorithm.}
  \item{unit.weights}{Constrains predictor loading weights to have unit norm.}
  \item{weight}{When \code{weight = TRUE}, a set of geometrically decaying weights is applied to the decentralized data.}
  \item{beta}{When \code{weight = TRUE}, the speed of geometric decay is governed by a weight parameter.}
  \item{interval}{When \code{interval = TRUE}, produces distributional forecasts.}
  \item{method}{Method used for computing prediction intervals.}
  \item{alpha}{\code{1-alpha} gives the nominal coverage probability.}
  \item{B}{Number of replications.}
  \item{adjust}{When \code{adjust = TRUE}, an adjustment is performed.}
  \item{backh}{When \code{adjust = TRUE}, an adjustment is performed by evaluating the difference between
               predicted and actual values in a testing set. \code{backh} specifies the testing set.}
}
\details{
\emph{\bold{Point forecasts:}}\cr

The NIPALS function implements the orthogonal scores algorithm, as described in Martens and Naes (1989).
This is one of the two classical PLSR algorthms, the other is the simple partial least squares regression in DeJong (1993).
The difference between these two approaches is that the NIPALS deflates the original predictors and responses,
while the SIMPLS deflates the covariance matrix of original predictors and responses. Thus, SIMPLS is more computationally efficient than NIPALS.

In a functional data set, the functional PLSR can be performed by setting the functional responses to be 1 lag ahead of the functional predictors.
This idea has been adopted from the Autoregressive Hilbertian processes of order 1 (ARH(1)) of Bosq (2000).

\emph{\bold{Distributional forecasts:}}\cr

\emph{Parametric method:}

Influenced by the works of Denham (1997) and Phatak et al. (1993), one way of constructing prediction intervals in the PLSR
is via a local linearization method (also known as the Delta method). It can be easily understood as the first two terms in a Taylor series
expansion. The variance of coefficient estimators can be approximated, from which an analytic-formula based prediction intervals are constructed.

\emph{Nonparametric method:}

After discretizing and decentralizing functional data \eqn{f_t(x)} and \eqn{g_s(y)}, a PLSR model with \eqn{K} latent components is built.
Then, the fit residuals \eqn{o_s(y_i)} between \eqn{g_s(y_i)} and \eqn{\hat{g}_s(y_i)}{hat(g)_s(y_i)} are calculated as
\deqn{o_s(y_i)=g_s(y_i)-\hat{g}_s(y_i), i=1,...,p.}{o_s(y_i)=[g_s(y_i)-hat(g)_s(y_i)], i=1,...,p.}

The next step is to generate \eqn{B} bootstrap samples \eqn{o_s^b(y_i)} by randomly sampling with replacement
from \eqn{[o_1(y_i),...,o_n(y_i)]}. Adding bootstrapped residuals to the original
response variables in order to generate new bootstrap responses,
\deqn{g_s^b(y_i)=g_s(y_i)+o_s^b(y_i).}\cr
Then, the PLSR models are constructed using the centered and discretized predictors and bootstrapped responses
to obtain the boostrapped regression coefficients and point forecasts, from which the empirical prediction intervals and kernel density plots are constructed.
}
\value{
A list containing the following components is returned.
  \item{B}{\eqn{(p \times m)}{(p x m)} matrix containing the regression coefficients. \eqn{p} is the number of variables in
  the predictors and \eqn{m} is the number of variables in the responses.}
  \item{P}{\eqn{(p \times order)}{(p x order)} matrix containing the predictor loadings.}
  \item{Q}{\eqn{(m \times order)}{(m x order)} matrix containing the response loadings.}
  \item{T}{\code{(ncol(data$y)-1) x order} matrix containing the predictor scores.}
  \item{R}{\eqn{(p\times order)}{(p x order)} matrix containing the weights used to construct the latent components of predictors.}
  \item{Yscores}{\code{(ncol(data$y)-1) x order} matrix containing the response scores.}
  \item{projection}{\eqn{(p\times order)}{(p x order)} projection matrix used to convert predictors to predictor scores.}
  \item{meanX}{An object of class \code{fts} containing the column means of predictors.}
  \item{meanY}{An object of class \code{fts} containing the column means of responses.}
  \item{Ypred}{An object of class \code{fts} containing the 1-step-ahead predicted values of the responses.}
  \item{fitted}{An object of class \code{fts} containing the fitted values.}
  \item{residuals}{An object of class \code{fts} containing the regression residuals.}
  \item{Xvar}{A vector with the amount of predictor variance explained by each number of component.}
  \item{Xtotvar}{Total variance in predictors.}
  \item{weight}{When \code{weight = TRUE}, a set of geometrically decaying weights is given. When \code{weight = FALSE}, weights are all equal 1.}
  \item{x1}{Time period of a \code{fts} object, which can be obtained from \code{colnames(data$y)}.}
  \item{y1}{Variables of a \code{fts} object, which can be obtained from \code{data$x}.}
  \item{ypred}{Returns the original functional predictors.}
  \item{y}{Returns the original functional responses.}
  \item{bootsamp}{Bootstrapped point forecasts.}
  \item{lb}{Lower bound of prediction intervals.}
  \item{ub}{Upper bound of prediction intervals.}
  \item{lbadj}{Adjusted lower bound of prediction intervals.}
  \item{ubadj}{Adjusted upper bound of prediction intervals.}
  \item{lbadjfactor}{Adjusted lower bound factor, which lies generally between 0.9 and 1.1.}
  \item{ubadjfactor}{Adjusted upper bound factor, which lies generally between 0.9 and 1.1.}
}
\references{
S. Wold and A. Ruhe and H. Wold and W. J. Dunn (1984) "The collinearity problem in linear regression. The partial least squares (PLS) approach to generalized inverses", \emph{SIAM Journal of Scientific and Statistical Computing}, \bold{5}(3), 735-743.

S. de Jong (1993) "SIMPLS: an alternative approach to partial least square regression", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{18}(3), 251-263.

C J. F. Ter Braak and S. de Jong (1993) "The objective function of partial least squares regression", \emph{Journal of Chemometrics}, \bold{12}(1), 41-54.

B. Dayal and J. MacGregor (1997) "Recursive exponentially weighted PLS and its applications to adaptive	control and prediction", \emph{Journal of Process Control}, \bold{7}(3), 169-179.

B. D. Marx (1996) "Iteratively reweighted partial least squares estimation for generalized linear regression", \emph{Technometrics}, \bold{38}(4), 374-381.

L. Xu and J-H. Jiang and W-Q. Lin and Y-P. Zhou and H-L. Wu and G-L. Shen and R-Q. Yu (2007) "Optimized sample-weighted partial least squares", \emph{Talanta}, \bold{71}(2), 561-566.

A. Phatak and P. Reilly and A. Penlidis (1993) "An approach to interval estimation in partial least squares regression", \emph{Analytica Chimica Acta}, \bold{277}(2), 495-501.

M. Denham (1997) "Prediction intervals in partial least squares", \emph{Journal of Chemometrics}, \bold{11}(1), 39-52.

D. Bosq (2000) \emph{Linear Processes in Function Spaces}, New York: Springer.

N. Faber (2002) "Uncertainty estimation for multivariate regression coefficients", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{64}(2), 169-179.

J. A. Fernandez Pierna and L. Jin and F. Wahl and N. M. Faber and D. L. Massart (2003) "Estimation of partial least squares regression prediction uncertainty when the reference values carry a sizeable measurement error", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{65}(2), 281-291.

P. T. Reiss and R. T. Ogden (2007), "Functional principal component regression and functional partial least squares", \emph{Journal of the American Statistical Association}, \bold{102}(479), 984-996.

C. Preda, G. Saporta (2005) "PLS regression on a stochastic process", \emph{Computational Statistics and Data Analysis}, \bold{48}(1), 149-158.

C. Preda, G. Saporta, C. Leveder (2007) "PLS classification of functional data", \emph{Computational Statistics}, \bold{22}, 223-235.

A. Delaigle and P. Hall (2012), "Methodology and theory for partial least squares applied to functional data", \emph{Annals of Statistics}, \bold{40}(1), 322-352.

M. Febrero-Bande, P. Galeano, W. Gonz\'alez-Manteiga (2017), "Functional principal component regression and functional partial least-squares regression: An overview and a comparative study", \emph{International Statistical Review}, \bold{85}(1), 61-83.
}
\author{Han Lin Shang}
\seealso{\code{\link[ftsa]{ftsm}}, \code{\link[ftsa]{forecast.ftsm}}, \code{\link[ftsa]{plot.fm}},
\code{\link[ftsa]{summary.fm}}, \code{\link[ftsa]{residuals.fm}}, \code{\link[ftsa]{plot.fmres}}}
\examples{
# When weight = FALSE, all observations are assigned equally.
# When weight = TRUE, all observations are assigned geometrically decaying weights.
fplsr(data = ElNino_ERSST_region_1and2, order = 6, type = "nipals")
fplsr(data = ElNino_ERSST_region_1and2, order = 6)
fplsr(data = ElNino_ERSST_region_1and2, weight = TRUE)
fplsr(data = ElNino_ERSST_region_1and2, unit.weights = FALSE)
fplsr(data = ElNino_ERSST_region_1and2, unit.weights = FALSE, weight = TRUE)

# The prediction intervals are calculated numerically.
fplsr(data = ElNino_ERSST_region_1and2, interval = TRUE, method = "delta")

# The prediction intervals are calculated by bootstrap method.
fplsr(data = ElNino_ERSST_region_1and2, interval = TRUE, method = "boota")
}
\keyword{models}
