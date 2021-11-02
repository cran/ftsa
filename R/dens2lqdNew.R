#' Function for converting densities to log quantile density functions
#'
#' @param dens density values on dSup - must be strictly positive (otherwise will truncate) and integrate to 1
#' @param dSup support (grid) for Density domain
#' @param lqdSup support of length M for lqd domain - must begin at 0 and end at 1; (default: seq(0, 1, length(dSup)))
#' @param t0 value in dSup for which the cdf value c is retained, i.e. c = F(t0) (default: dSup[1])
#' @param verbose if FALSE, repress some messages (default: TRUE)
#'
#' @return list with 'lqdSup' a grid on [0,1], 'lqd' the log quantile density on lqdSup, and 'c' the value of the cdf at t0
#'
#' @examples
#' x <- seq(0,2,length.out =512)
#' y <- rep(0.5,length.out =512)
#' y.lqd <- dens2lqd( dens=y, dSup = x) # should equate # log(2)
#' @references
#' \cite{Functional Data Analysis for Density Functions by Transformation to a Hilbert space, Alexander Petersen and Hans-Georg Mueller, 2016}
#' @export

dens2lqd = function(dens, dSup, lqdSup = seq(0, 1, length.out = length(dSup)), t0 = dSup[1], verbose = TRUE){

  # Check density requirements
  if(any(dens < 0)){
    stop('Please correct negative density values.')
  }

  if(abs( trapzRcpp(X = dSup, Y = dens) - 1) > 1e-5){
    
    warning('Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.')
    dens = dens/trapzRcpp(X = dSup, Y = dens)
    
  }
  
  if(any(dens == 0)){
    if(verbose){
    print("There are some zero density values - truncating support grid so all are positive")
    }
    lbd = min(which(dens > 0))
    ubd = max(which(dens > 0))
    dens = dens[lbd:ubd]
    dSup = dSup[lbd:ubd]
    dens = dens/trapzRcpp(X = dSup, Y = dens)
  }
  
  N = length(dSup)
  
  # Check LQD output grid
  if(is.null(lqdSup)){
      
    lqdSup = seq(0, 1, length.out = N)
  
  }else if(!all.equal( range(lqdSup),c(0,1) )){
    
    if(verbose){
    print("Problem with support of the LQD domain's boundaries - resetting to default.")
    }
    lqdSup = seq(0, 1, length.out = N)

  }
  
  # Check t0
  if(!(t0 %in% dSup)){
    
    if(verbose){
      print("t0 is not a value in dSup - resetting to closest value")
    }
    t0 = dSup[which.min(abs(dSup - t0))]
    
  }
  
  M = length(lqdSup) 
  c_ind = which(dSup == t0)

  # Get CDF and lqd on temporary grid, compute c
  tmp = cumtrapzRcpp(X = dSup, dens)
  c = tmp[c_ind]

  indL = duplicated(tmp[1:floor(N/2)])
  indR = duplicated(tmp[(floor(N/2)+1):N], fromLast = TRUE)
  qtemp = tmp[!c(indL, indR)]
  lqd_temp = -log(dens[!c(indL, indR)]);
  
  # Interpolate lqdSup, keeping track of Inf values at boundary, then compute c
  lqd = rep(0, 1, M)

  if(any(is.infinite(lqd_temp[c(1, N)]))){

    tmpInd = 1:N
    Ind = 1:M
    
    if(lqd_temp[1] == Inf){
      
      lqd[1] = Inf
      tmpInd = tmpInd[-1]
      Ind = Ind[-1]
    
    }
    
    if(lqd_temp[N] == Inf){
  
      lqd[M] = Inf
      tmpInd = tmpInd[-length(tmpInd)]
      Ind = Ind[-length(Ind)]
      
    }
    
    lqd[Ind] = approx(x = qtemp[tmpInd], y = lqd_temp[tmpInd], xout = lqdSup[Ind], rule = 2)$y 
    
  }else{
    
    lqd = approx(x = qtemp, y = lqd_temp, xout = lqdSup, rule = c(2,2))$y 
    
  }
  
  return(list('lqdSup',  lqdSup, 'lqd' = lqd, 'c' = c))
}
