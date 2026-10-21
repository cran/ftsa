

####################################
### Fast fractional differencing ###
####################################

#' @title Fast fractional differencing procedure of Jensen and Nielsen (2014).
#' @description Takes the d-th difference of the series.
#' @details This code was first published on the \href{http://www.econ.queensu.ca/faculty/mon/software/}{university webpage of Morten Nielsen} and is 
#' redistributed here with the author's permission.
#' @param x series to be differenced
#' @param d memory parameter indicating order of the fractional difference.
#' @author Jensen, A. N. and Nielsen, M. O.
#' @examples acf(fdiff(x=rnorm(500), d=0.4))
#' @references Jensen, A. N. and Nielsen, M. O. (2014):
#' A fast fractional difference algorithm, Journal of Time Series Analysis 35(5), pp. 428-436.
#' @export

fdiff<- function(x, d){
  iT <- length(x)
  np2 <- nextn(2*iT - 1, 2)
  k <- 1:(iT-1)
  b <- c(1, cumprod((k - d - 1)/k))
  dx <- fft(fft(c(b, rep(0, np2-iT)))*
              fft(c(x, rep(0, np2-iT))), inverse=T)/np2;
  return(Re(dx[1:iT]))
}