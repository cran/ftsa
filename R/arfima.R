## Automatic ARFIMA modelling
## Will return Arima object if d < 0.01 to prevent estimation problems
arfima <- function(x, drange = c(0, 0.5), ...)
{
    # Remove mean
    meanx <- mean(x)
    x <- x - meanx
    
    # Choose differencing parameter with AR(2) proxy to handle correlations
    warn <- options(warn=-1)$warn
    fit <- fracdiff(x,nar=2)
    options(warn=warn)
   
    # Choose p and q
    d <- fit$d
    y <- diffseries(x, d=d)
    fit <- auto.arima(y, max.P=0, max.Q=0, stationary=TRUE, ...)
    
    # Refit model using fracdiff
    warn <- options(warn=-1)$warn
    fit <- fracdiff(x, nar=fit$arma[1], nma=fit$arma[2])
    options(warn=warn)
    
    # Add things to model that will be needed by forecast.fracdiff
    fit$x <- x + meanx
    fit$residuals <- residuals(fit)
    fit$fitted <- fit$x - fit$residuals
    
    return(fit)
}

