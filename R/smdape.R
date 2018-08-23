smdape <- function(forecast, true)
{
     if (length(forecast) != length(true))
         stop("SMDAPE: the lengths of input vectors must be the same.")
     err = median(200 * abs(true - forecast) / (true + forecast))
     return(round(err, 6))
}
