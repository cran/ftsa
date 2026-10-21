
#---------------------------       profiled/ concentrated likelihood function       ------------------------------------#

#' Concentrated local Whittle likelihood. Only for internal use. cf. Shimotsu and Phillips (2005), p. ???.
#'@keywords internal
wd.elw<-function(d){1/2*(1+cos(4*pi*d))}

#' concentrated likelihood function for ELW estimator
#'@keywords internal
R.elw<-function(d,data,m){
  T<-length(data)
  lambda<-2*pi/T
  Peri<-per(fdiff(data, d=d))[-1]
  K<-log(mean(Peri[1:m]))-2*d*mean(log(lambda*(1:m)))
  K
}

#' concentrated likelihood function for ELW estimator - weighted version
#'@keywords internal
R.elw.weighted<-function(d,data,m){
  data<-(data-wd.elw(d)*mean(data)-(1-wd.elw(d))*data[1])[-1]
  T<-length(data)
  lambda<-2*pi/T
  Peri<-per(fdiff(data, d=d))[-1]
  K<-log(mean(Peri[1:m]))-2*d*mean(log(lambda*(1:m)))
  K
}

#---------------------------      Two-Step Exact local whittle estimation       ------------------------------------#

#' @title Two-Step Exact local Whittle estimator of fractional integration with 
#' unknown mean and time trend.
#' @description \code{ELW2S} implements the two-step ELW estimator of 
#' Shimotsu (2010) that is consistent and asymptotically normal in the range from -1/2 to 2.
# #' @details add details here
#' @param data data vector of length T.
#' @param m bandwith parameter specifying the number of Fourier frequencies.
#' used for the estimation usually \code{floor(1+T^delta)}, where 0<delta<1.
#' @param trend_order specifies the form of detrending: 0 for a constant, only, 
#' 1 for a linear trend, and so on. 
#' @param taper string from \code{c("Velasco","HC")} specifying the tapered form of the local Whittle estimator used in the first step.
#' @author Christian Leschinski
#' @references Shimotsu, K. (2010): Exact Local Whittle
#' Estimation Of Fractional Integration with Unknown Mean and Time Trend. Econometric Theory,
#'  Vol. 26, pp. 501 - 540.
#' @examples
#' library(fracdiff)
#' T<-1000
#' d<-0.8
#' trend<-(1:T)/T
#' series<-cumsum(fracdiff.sim(T,d=(d-1))$series)
#' ts.plot(series)
#' ELW2S(series, m=floor(1+T^0.7), trend_order=0)$d
#' series2<-series+2*trend
#' ELW2S(series2, m=floor(1+T^0.7), trend_order=1)$d
#' series3<-series+2*trend+2*trend^2
#' ELW2S(series3, m=floor(1+T^0.7), trend_order=2)$d
#' @export
ELW2S<-function(data, m, trend_order=0, taper=c("Velasco","HC")){
  taper<-taper[1]
  aux_est<-local.W(data=data, m=m, taper=taper, int=c(-1/2,2.5))
  d_init<-aux_est$d
  T<-length(data)
  if(trend_order==0){
    Xt<-residuals(lm(data~1))}else{
      time<-1:T
      Xt<-residuals(lm(data~poly(time,trend_order)))
    }
  d.hat<-optim(par=d_init, fn=R.elw.weighted, method="BFGS", data=data, m=m)$par
  if(abs(d.hat-d_init)>2.5*aux_est$s.e.){d.hat<-d_init}
  se<-1/(2*sqrt(m))
  return(list("d"=d.hat, "s.e."=se))
}

