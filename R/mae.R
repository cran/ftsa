mae = function(forecast, true)
{
    if (length(forecast) != length(true))
        stop("MAE: the lengths of input vectors must be the same.")
    err = mean(abs(forecast - true))
    return(round(err, 6))
}

