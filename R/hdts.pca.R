## high dimensional time series pca

hdpca <- function(x, order, h0 = 5)
{
  xs <- scale(x, scale = FALSE)
  n <- nrow(x) # n is sample size
  m <- ncol(x) # m is number of populations
  c1 <- array(0, dim = c(n-1, m, m)) # c1[h,,] is the covariance at h
  for(h in 1: (n-1))
  {
      for(k in 1: (n-h))
      {
          c1[h,,] <- c1[h,,] + xs[k+h,]%*%t(xs[k,])
      }
      c1[h,,] <- c1[h,,]/(n-h)
  }
  # choose q 
  f <- array(0, dim = c(m, m)) # f is weighted sum of the covariance at all lags
  for (h in 1:h0)
  {
      f <- f + c1[h,,]%*%t(c1[h,,])
  }
  md <- eigen(f, symmetric = TRUE)
  md$values[md$values<0] <- 0
  score <- xs%*%(md$vectors[, 1:order])
  fitted <- t(score%*%t(md$vectors[, 1:order])) + apply(x, 2, mean)
  varprop <- cumsum(md$values/sum(md$values))[1:order]
  basis <- cbind(apply(x, 2, mean), md$vectors[, 1:order])
  resid <- x-t(fitted)
  return(list(coef= score, fitted =fitted, basis = basis, varprop =varprop, order = order))
}


  
