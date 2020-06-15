\name{sim_ex_cluster}
\alias{sim_ex_cluster}
\alias{sim_ex_cluster.smooth}
\docType{data}
\title{
Simulated multiple sets of functional time series
}
\description{
We generate  2 groups of \eqn{m} functional time series. For each \eqn{i} in \{1, \dots, m\} in a given cluster \eqn{c}, \eqn{c} in \{1,2\}, the \eqn{t} th function, \eqn{t} in \{1,\dots, T\}, is given by
\deqn{Y_{it}^{(c)} (x)= \mu^{(c)}(x) + \sum_{k=1}^{2}\xi_{tk}^{(c)} \rho_k^{(c)} (x) + \sum_{l=1}^{2}\zeta_{itl}^{(c)} \psi_l^{(c)} (x) + \upsilon_{it}^{(c)} (x)}
}
\usage{data("sim_ex_cluster")}
\details{
The mean functions for each of these two clusters are set to be \eqn{\mu^{(1)}(x) = 2(x-0.25)^{2}} and \eqn{\mu^{(2)}(x) = 2(x-0.4)^{2}+0.1}. 

While the variates \eqn{\mathbf{\xi_{tk}^{(c)}}=(\xi_{1k}^{(c)}, \xi_{2k}^{(c)}, \ldots, \xi_{Tk}^{(c)})^{\top}} for both clusters, are generated from autoregressive of order 1 with parameter 0.7, while the variates \eqn{\zeta_{it1}^{(c)}} and \eqn{\zeta_{it2}^{(c)}} for both clusters, are generated from independent and identically distributed \eqn{N(0,0.5)} and \eqn{N(0,0.25)}, respectively.

The basis functions for the common-time trend for the first cluster, \eqn{\rho_k^{(1)} (x)}, for \eqn{k} in \{1,2\} are \eqn{sqrt(2)*sin(\pi*(0:200/200))} and \eqn{sqrt(2)*cos(\pi*(0:200/200))} respectively; and the basis functions for the common-time trend for the second cluster, \eqn{\rho_k^{(2)} (x)}, for \eqn{k} in \{1,2\} are \eqn{sqrt(2)*sin(2\pi*(0:200/200))} and \eqn{sqrt(2)*cos(2\pi*(0:200/200))} respectively.

The basis functions for the residual for the first cluster, \eqn{\psi_l^{(1)} (x)}, for \eqn{l} in \{1,2\} are \eqn{sqrt(2)*sin(3\pi*(0:200/200))} and \eqn{sqrt(2)*cos(3\pi*(0:200/200))} respectively; and the basis functions for the residual for the second cluster, \eqn{\psi_l^{(2)} (x)}, for \eqn{l} in \{1,2\} are \eqn{sqrt(2)*sin(4\pi*(0:200/200))} and \eqn{sqrt(2)*cos(4\pi*(0:200/200))} respectively.

The measurement error \eqn{\upsilon_{it}} for each continuum x is generated from independent and identically distributed \eqn{N(0, 0.2^2)}
}
\examples{
data(sim_ex_cluster)
}
\keyword{datasets}
