sco.resamp <- function(object, h, B)
{
    n <- length(object$x)
    sco <- object$x
    olivia <- object$mean
    
    forerr = list()
    for(i in 1: h)
    {
        forerr[[i]] = array(NA, (n/2-h+1))
    }
    for(j in 1:h)
    {
        for(i in 1:(n/2-j+1))
        {
            fore = forecast(auto.arima(sco[1:(n/2-1+i)]), h = j)$mean[j]
            forerr[[j]][i] = sco[n/2-1+i+j] - fore
        }
    }
    
    ny = array(NA, dim = c(B, h))
    for(j in 1:h)
    {
        ny[ , j] = sample(forerr[[j]], size = B, replace = TRUE)
    }
    oli <- t(array(rep(olivia, B * h), dim = c(h, B)))
    fo <- oli + ny
    return(fo)
}
