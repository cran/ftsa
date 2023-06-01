forecast.dfpca <- function(object, h = 10, ...)
{
    k <- object$order
    y <- array(NA, dim = c(h, k))
    for(ik in 1:k)
    {
        y[, ik] <- forecast(auto.arima(object$coef[, ik]), h = h)$mean
    }
    f <- object$basis[, 2: (k+1)]%*%t(y) + object$basis[,1]
    return(f)
}
