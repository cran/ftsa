dfpca <- function(x, order, q)
{
    xs <- scale(x, scale = FALSE)
    n <- nrow(x) # n is sample size
    m <- ncol(x) # m is age span
    c1 <- array(0, dim = c(n, m, m)) # c1[h,,] is the covariance at h-1
    c2 <- array(0, dim = c(n, m, m)) # c2[h,,] is the covariance at 1-h
    for(h in 1: n)
    {
        for(k in h:n)
        {
            c1[h,,] <- c1[h,,] + xs[k,]%*%t(xs[k-h+1,])
            c2[h,,] <- c2[h,,] + xs[k-h+1,]%*%t(xs[k,])
        }
        c1[h,,] <- c1[h,,]/n
        c2[h,,] <- c2[h,,]/n
    }
    # choose q
    f <- array(0, dim = c(m, m)) # f is weighted sum of the covariance at all lags
    for(h in 2:n)
    {
        #  f <- f + wt.tri(h/q)*c1[h,,] + wt.tri(-h/q)*c2[h,,]
        f <- f + (wtflat((h/q), 0.5))*c1[h,,] + (wtflat((-h/q), 0.5))*c2[h,,]
    }
    f <- f + c1[1,,]
  
    md <- eigen(f, symmetric = TRUE)
    md$values[md$values<0] <- 0
    score <- xs%*%md$vectors[, 1:order]
    fitted <- t(score%*%t(md$vectors[, 1:order])) + apply(x, 2, mean)
    varprop <- cumsum(md$values/sum(md$values))[1:order]
    basis <- cbind(apply(x, 2, mean), md$vectors[, 1:order])
    return(list(coef= score, fitted =fitted, basis = basis, varprop =varprop, order = order))
}

# weight function (flat top)

wtflat <- function(x, c)
{
    if(-1< x & x <= -c)
    {
        return(x/(1-c) + 1/(1-c))
    }
    if(-c< x & x < c)
    {
        return(1)
    }
    if( c< x & x <1)
    {
        return(x/(c-1) - 1/(c-1))
    }
    else return(0)
}

# weight function (triangle)

wt.tri <- function(x)
{
    if(abs(x)>1)
    {
        return(0)
    }
    else
    {
        return(1-abs(x))
    }
}

forecast.dfpca <- function(object, h = 10)
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
