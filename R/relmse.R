relmse = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("RelMAE: the lengths of input vectors must be the same.")
  ferror = ftsa:::mse(forecast, true)
  ferrorbench = ftsa:::mse(forecastbench, true)
  relativerror = ferror / ferrorbench
  return(round(relativerror, 6))
}
