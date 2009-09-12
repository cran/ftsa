gmrae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("GMRAE: the lengths of input vectors must be the same.")
  ferror = (forecast - true)
  ferrorbench = (forecastbench - true)
  relativerror = prod(abs(ferror/ferrorbench))^(1/length(forecast))
  return(round(relativerror, 6))
}

