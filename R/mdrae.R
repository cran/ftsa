mdrae = function(forecast, forecastbench, true)
{
  if (length(forecast) != length(true))
      stop("MdRAE: the lengths of input vectors must be the same.")
  ferror = (forecast - true)
  ferrorbench = (forecastbench - true)
  relativerror = median(abs(ferror / ferrorbench))
  return(round(relativerror, 6))
}
  