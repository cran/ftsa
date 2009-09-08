mase = function(forecast,true)
{ 
  if (length(forecast) != length(true))
      stop("MASE: the lengths of input vectors must be the same.")
  n = length(true)
  truerror = vector(,n-1)
  for(i in 1:(n-1)){
      truerror[i] = abs(true[i+1] - true[i])
  }
  qt = (forecast-true) / (sum(truerror) / (n-1))
  scalederror = mean(abs(qt))
  return(round(scalederror, 6))
}