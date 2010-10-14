relmae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("RelMAE: the lengths of input vectors must be the same.")
  ferror = ftsa:::mae(forecast, true)
  ferrorbench = ftsa:::mae(forecastbench, true)
  relativerror = ferror / ferrorbench
  return(round(relativerror, 6))
}
