\name{ftsa-package}
\alias{ftsa-package}
\alias{ftsa}
\docType{package}
\title{
Functional Time Series Analysis
}
\description{
This package presents descriptive statistics of functional data; implements principal component regression and partial least squares regression to provide point and distributional forecasts for functional data; utilizes functional linear regression, ordinary least squares, penalized least squares, ridge regression, and moving block approaches to dynamically update point and distributional forecasts when partial data points in the most recent curve are observed; performs stationarity test for a functional time series;
    estimates a long-run covariance function by kernel sandwich estimator.
}
\author{
Rob J Hyndman and Han Lin Shang 

Maintainer: Han Lin Shang <hanlin.shang@anu.edu.au>
}
\references{
###########################
# References in Statistics
###########################

R. J. Hyndman and H. L. Shang (2009) "Forecasting functional time series (with discussion)", \emph{Journal of the Korean Statistical Society}, \bold{38}(3), 199-221.

R. J. Hyndman and H. L. Shang (2010) "Rainbow plots, bagplots, and boxplots for functional data", \emph{Journal of Computational and Graphical Statistics}, \bold{19}(1), 29-45.

H. L. Shang and R. J. Hyndman (2011) "Nonparametric time series forecasting with dynamic updating", \emph{Mathematics and Computers in Simulation}, \bold{81}(7), 1310-1324.

H. L. Shang (2011) "rainbow: an R package for visualizing functional time series, \emph{The R Journal}, \bold{3}(2), 54-59.

H. L. Shang (2013) "Functional time series approach for forecasting very short-term electricity demand", \emph{Journal of Applied Statistics}, \bold{40}(1), 152-168.

H. L. Shang (2013) "ftsa: An R package for analyzing functional time series", \emph{The R Journal}, \bold{5}(1), 64-72.

H. L. Shang (2014) "A survey of functional principal component analysis", \emph{Advances in Statistical Analysis}, \bold{98}(2), 121-142.

H. L. Shang (2014) "Bayesian bandwidth estimation for a functional nonparametric regression model with mixed types of regressors and unknown error density", \emph{Journal of Nonparametric Statistics}, \bold{26}(3), 599-615.

H. L. Shang (2014) "Bayesian bandwidth estimation for a semi-functional partial linear regression model with unknown error density", \emph{Computational Statistics}, \bold{29}(3-4), 829-848.

H. L. Shang (2015) "Resampling techniques for estimating the distribution of descriptive statistics of functional data", \emph{Communications in Statistics - Simulation and Computation}, \bold{44}(3), 614- 635.

H. L. Shang (2016) "Mortality and life expectancy forecasting for a group of populations in developed countries: A robust multilevel functional data method", in C. Agostinelli, A. Basu, P. Filzmoser, D. Mukherjee (ed.), Recent Advances in Robust Statistics: Theory and Applications, Springer, India, pp. 169-184.

H. L. Shang (2016) "Mortality and life expectancy forecasting for a group of populations in developed countries: A multilevel functional data method", \emph{Annals of Applied Statistics}, \bold{10}(3), 1639-1672. 

H. L. Shang (2016) "A Bayesian approach for determining the optimal semi-metric and bandwidth in scalar-on-function quantile regression with unknown error density and dependent functional data", \emph{Journal of Multivariate Analysis}, \bold{146}, 95-104.

H. L. Shang (2017) "Functional time series forecasting with dynamic updating: An application to intraday particulate matter concentration", \emph{Econometrics and Statistics}, \bold{1}, 184-200.

H. L. Shang (2017) "Forecasting Intraday S&P 500 Index Returns: A Functional Time Series Approach", \emph{Journal of Forecasting}, \bold{36}(7), 741-755.

H. L. Shang and R. J. Hyndman (2017) "Grouped functional time series forecasting: An application to age-specific mortality rates", \emph{Journal of Computational and Graphical Statistics}, \bold{26}(2), 330-343.

G. Rice and H. L. Shang (2017) "A plug-in bandwidth selection procedure for long-run covariance estimation with stationary functional time series", \emph{Journal of Time Series Analysis}, \bold{38}(4), 591-609.

P. Reiss, J. Goldsmith, H. L. Shang and R. T. Ogden (2017) "Methods for scalar-on-function regression", \emph{International Statistical Review}, \bold{85}(2), 228-249.

P. Kokoszka, G. Rice and H. L. Shang (2017) "Inference for the autocovariance of a functional time series under conditional heteroscedasticity", \emph{Journal of Multivariate Analysis}, \bold{162}, 32-50.

Y. Gao, H. L. Shang and Y. Yang (2017) "High-dimensional functional time series forecasting", in G. Aneiros, E. Bongiorno, R. Cao and P. Vieu (ed.), Functional Statistics and Related Fields, Springer, Cham, pp. 131-136.

Y. Gao and H. L. Shang (2017) "Multivariate functional time series forecasting: An application to age-specific mortality rates", \emph{Risks}, \bold{5}(2), Article 21.

H. L. Shang (2018) "Visualizing rate of change: An application to age-specific fertility rates", \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, \bold{182}(1), 249-262.

H. L. Shang (2018) "Bootstrap methods for stationary functional time series", \emph{Statistics and Computing}, \bold{28}(1), 1-10.

Y. Gao, H. L. Shang and Y. Yang (2019) "High-dimensional functional time series forecasting: An application to age-specific mortality rates", \emph{Journal of Multivariate Analysis}, \bold{170}, 232-243.

D. Li, P. M. Robinson and H. L. Shang (2020) "Long-range dependent curve time series", \emph{Journal of the American Statistical Association: Theory and Methods}, \bold{115}(530), 957-971.

H. L. Shang (2020) "A comparison of Hurst exponent estimators in long-range dependent curve time series", \emph{Journal of Time Series Econometrics}, \bold{12}(1).

D. Li, P. M. Robinson and H. L. Shang (2021) "Local Whittle estimation of long range dependence for functional time series", \emph{Journal of Time Series Analysis}, \bold{42}(5-6), 685-695.

H. L. Shang and R. Xu (2021) "Functional time series forecasting of extreme values", \emph{Communications in Statistics: Case Studies, Data Analysis and Applications}, \bold{7}(2), 182-199.

U. Beyaztas, H. L. Shang and Z. Yaseen (2021) "Development of functional autoregressive model based exogenous hydrometeorological variables for river flow prediction", \emph{Journal of Hydrology}, \bold{598}, 126380.

U. Beyaztas and H. L. Shang (2022) "Machine learning-based functional time series forecasting: Application to age-specific mortality rates", \emph{Forecasting}, \bold{4}(1), 394-408.

Y. Yang, Y. Yang and H. L. Shang (2022) "Feature extraction for functional time series: Theory and application to NIR spectroscopy data", \emph{Journal of Multivariate Analysis}, \bold{189}, 104863.

A. E. Fernandez, R. Jimenez and H. L. Shang (2022) "On projection methods for functional time series forecasting", \emph{Journal of Multivariate Analysis}, \bold{189}, 104890.

H. L. Shang (2022) "Not all long-memory estimators are born equal: A case of non-stationary curve time series", \emph{The Canadian Journal of Statistics}, \bold{50}(1), 357-380.

X. Huang, H. L. Shang and D. Pitt (2022) "Permutation entropy and its variants for measuring temporal dependence", \emph{Australian and New Zealand Journal of Statistics}, \bold{64}(4), 442-477.

H. L. Shang, J. Cao and P. Sang (2022) "Stopping time detection of wood panel compression: A functional time series approach", \emph{Journal of the Royal Statistical Society: Series C}, \bold{71}(5), 1205-1224.

C. Tang, H. L. Shang and Y. Yang (2022) "Clustering and forecasting multiple functional time series", \emph{The Annals of Applied Statistics}, \bold{16}(4), 2523-2553.

J. Trinka, H. Haghbin, M. Maadooliat and H. L. Shang (2023) "Functional time series forecasting: Functional singular spectrum analysis approaches", \emph{Stat}, \bold{12}(1), e621.

D. Li, P. M. Robinson and H. L. Shang (2023) "Nonstationary fractionally integrated functional time series", \emph{Bernoulli}, \bold{29}(2), 1505-1526.

X. Huang and H. L. Shang (2023) "Nonlinear autocorrelation function of functional time series", \emph{Nonlinear Dynamics: An International Journal of Nonlinear Dynamics and Chaos in Engineering Systems}, \bold{111}, 2537-2554.

H. L. Shang (2023) "Sieve bootstrapping memory parameter in long-range dependent stationary functional time series", \emph{AStA Advances in Statistical Analysis}, \bold{107}, 421-441.

E. Paparoditis and H. L. Shang (2023) "Bootstrap prediction bands for functional time series", \emph{Journal of the American Statistical Association: Theory and Methods}, \bold{118}(542), 972-986.

Y. Gao, H. L. Shang and Y. Yang (2024) "Factor-augmented smoothing model for functional data", \emph{Statistica Sinica}, \bold{34}(1), 1-26.

#############################
# References in Population Studies
#############################

H. L. Shang, H. Booth and R. J. Hyndman (2011) "Point and interval forecasts of mortality rates and life expectancy: a comparison of ten principal component methods, \emph{Demographic Research}, \bold{25}(5), 173-214.

H. L. Shang (2012) "Point and interval forecasts of age-specific fertility rates: a comparison of functional principal component methods", \emph{Journal of Population Research}, \bold{29}(3), 249-267.

H. L. Shang (2012) "Point and interval forecasts of age-specific life expectancies: a model averaging", \emph{Demographic Research}, \bold{27}, 593-644.

H. L. Shang, A. Wisniowski, J. Bijak, P. W. F. Smith and J. Raymer (2014) "Bayesian functional models for population forecasting", in M. Marsili and G. Capacci (eds), Proceedings of the Sixth Eurostat/UNECE Work Session on Demographic Projections, Istituto nazionale di statistica, Rome, pp. 313-325.

H. L. Shang (2015) "Selection of the optimal Box-Cox transformation parameter for modelling and forecasting age-specific fertility", \emph{Journal of Population Research}, \bold{32}(1), 69-79.

H. L. Shang (2015) "Forecast accuracy comparison of age-specific mortality and life expectancy: Statistical tests of the results", \emph{Population Studies}, 69(3), 317-335.

H. L. Shang, P. W. F. Smith, J. Bijak, A. Wisniowski (2016) "A multilevel functional data method for forecasting population, with an application to the United Kingdom, \emph{International Journal of Forecasting}, 32(3), 629-649.

H. L. Shang (2017) "Reconciling forecasts of infant mortality rates at national and sub-national levels: Grouped time-series method", \emph{Population Research and Policy Review}, 36(1), 55-84.

R. J. Hyndman, Y. Zeng and H. L. Shang (2021) "Forecasting the old-age dependency ratio to determine the best pension age", \emph{Australian and New Zealand Journal of Statistics}, \bold{63}(2), 241-256.

Y. Yang and H. L. Shang (2022) "Is the group structure important in grouped functional time series?", \emph{Journal of Data Science}, \bold{20}(3), 303-324.

H. L. Shang and Y. Yang (2022) "Forecasting Australian subnational age-specific mortality rates", \emph{Journal of Population Research}, \bold{38}, 1-24.

Y. Yang, H. L. Shang and J. Raymer (2024) "Forecasting Australian fertility by age, region, and birthplace", \emph{International Journal of Forecasting}, in press.

###########################
# References in Actuarial Studies
###########################

H. L. Shang and S. Haberman (2017) "Grouped multivariate and functional time series forecasting: An application to annuity pricing", Presented at the Living to 100 Symposium, Orlando Florida, January 4-6, 2017.

H. L. Shang and S. Haberman (2017) "Grouped multivariate and functional time series forecasting: An application to annuity pricing", \emph{Insurance: Mathematics and Economics}, \bold{75}, 166-179.

H. L. Shang and S. Haberman (2018) "Model confidence sets and forecast combination: An application to age-specific mortality", \emph{Genus - Journal of Population Sciences}, \bold{74}, Article number: 19.

H. L. Shang and S. Haberman (2020) "Forecasting multiple functional time series in a group structure: an application to mortality", \emph{ASTIN Bulletin}, \bold{50}(2), 357-379.

H. L. Shang (2020) "Dynamic principal component regression for forecasting functional time series in a group structure", \emph{Scandinavian Actuarial Journal}, \bold{2020}(4), 307-322.

H. L. Shang and S. Haberman (2020) "Forecasting age distribution of death counts: An application to annuity pricing", \emph{Annals of Actuarial Science}, \bold{14}(1), 150-169.

H. L .Shang and S. Haberman and R. Xu (2022) "Multi-population modelling and forecasting age-specific life-table death counts", \emph{Insurance: Mathematics and Economics}, \bold{106}, 239-253.

####################
# References in Finance
####################

F. Kearney and H. L. Shang (2020) "Uncovering predictability in the evolution of the WTI oil futures curve", \emph{European Financial Management}, \bold{26}(1), 238-257.

H. L. Shang, K. Ji and U. Beyaztas (2021) "Granger causality of bivariate stationary curve time series", \emph{Journal of Forecasting}, \bold{40}(4), 626-635.

S. Butler, P. Kokoszka, H. Miao and H. L. Shang (2021) "Neural network prediction of crude oil futures using B-splines", \emph{Energy Economics}, \bold{94}, 105080.

H. L. Shang and F. Kearney (2022) "Dynamic functional time series forecasts of foreign exchange implied volatility surfaces", \emph{International Journal of Forecasting}, \bold{38}(3), 1025-1049. 

H. L. Shang and K. Ji (2023) "Forecasting intraday financial time series with sieve bootstrapping and dynamic updating", \emph{Journal of Forecasting}, \bold{42}(8), 1973-1988.
}
\keyword{package}

