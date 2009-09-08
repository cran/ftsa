rmsse = function(forecast, true)
{ 
  if (length(forecast) != length(true))
      stop("RMSSE: the lengths of input vectors must be the same.")
  n = length(true)
  truerror = vector(, n-1)
  for(i in 1:(n-1)){
      truerror[i] = (true[i+1] - true[i])^2
  }
  qt = (forecast - true) / (sum(truerror) / (n-1))
  scalederror = (mean(qt))^(1/2)
  return(round(scalederror, 6))
}
