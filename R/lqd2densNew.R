#' Function for converting log quantile densities to densities
#' 
#' @param  lqd log quantile density on lqdSup
#' @param  lqdSup support for lqd domain - must begin at 0 and end at 1
#' @param  t0 value for which the target cdf has F(t0) = c (default: 0)
#' @param  c value in lqdSup representing the value of target cdf at t0 (default: lqdSup[1])
#' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
#' @param  cut vector with two elements, indicating how many boundary to cut on the left and right side (default: c(0, 0)).  More will be cut off if exp(lqd) is infinite for some values.
#' @param  verbose if FALSE, repress some messages (default: TRUE)
#' 
#' @return list with 'dSup' the support grid and 'dens' the density values on dSup.  Due to cutoffs, 'dens' may integrate to less than one.
#' 
#' @examples
#' x <- seq(0,2,length.out =512)
#' y.lqd <- rep(log(2), times = 512)
#' y <- lqd2dens(lqd = y.lqd) # should equate # 1/2
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016} 
#' @export

lqd2dens = function(lqd, lqdSup = seq(0, 1, length.out = length(lqd)), dSup, t0 = 0, c = 0, useSplines = TRUE, cut = c(0, 0), verbose = TRUE){

  if(!all.equal( range(lqdSup),c(0,1) )){

    warning("Problem with support of the LQD domain's boundaries - resetting to default.")
    lqdSup = seq(0, 1, length.out = length(lqd))
    
  }
  
  M = length(lqd)
  r = which(exp(lqd) == Inf)
  
  if(length(r) > 0){
    
    if(any(r < floor(M/2))){
    cut[1] = max(cut[1], max(r[r < floor(M/2)]))
    }
    if(any(r >= floor(M/2))){
      cut[2] = max(cut[2], M - min(r[r >= floor(M/2)]) + 1)
    }    
  }

# Cut boundaries
  lqdSup = lqdSup[(cut[1] + 1):(M - cut[2])]
  lqd = lqd[(cut[1] + 1):(M - cut[2])]
  M = length(lqd) # reset N
  
  if(!(c %in% lqdSup)){
  
    if(c < lqdSup[1] || c > lqdSup[M]){
      
      stop("c is not contained withing range of lqdSup after cutoff")
    
    }
    
    if(verbose){
      
      print("c is not equal to a value in lqdSup - resetting to closest value")
    
    }
    c = lqdSup[which.min(abs(lqdSup - c))]
  
  }

  c_ind = which(lqdSup == c)
  
  if( useSplines ){    # Could fit spline if this yields more accurate numerical integration

    lqd_sp = splinefun(lqdSup, lqd, method = 'natural')
    lqd_exp = function(t) exp(lqd_sp(t))
    # Get grid for density space
    dtemp = t0 + c(0, cumsum(sapply(2:length(lqdSup), function(i) integrate(lqd_exp, lqdSup[i - 1], lqdSup[i])$value))) - integrate(lqd_exp, lqdSup[1], lqdSup[c_ind])$value

  } else {
    # Get grid and function for density space
    dtemp = t0 + cumtrapzRcpp(lqdSup, exp(lqd)) - trapzRcpp(lqdSup[1:c_ind], exp(lqd[1:c_ind]))
  }

  # Remove duplicates
  indL = duplicated(dtemp[1:floor(M/2)], fromLast = TRUE)
  indR = duplicated(dtemp[(floor(M/2)+1):M])
  dtemp = dtemp[!c(indL, indR)]
  dens_temp = exp(-lqd[!c(indL, indR)]);
  
  # Interpolate to dSup and normalize
  dSup = seq(dtemp[1], dtemp[length(dtemp)], length.out = M)
  dens = approx(x = dtemp, y = dens_temp, xout = dSup, rule = c(2,2))[[2]]
  dens = dens/trapzRcpp(X = dSup,Y = dens)*(lqdSup[M] - lqdSup[1]); # Normalize, accounting for boundary cutoff

  return(list('dSup' = dSup, 'dens' = dens))
}
