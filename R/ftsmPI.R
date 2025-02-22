ftsmPI <- function (object, B, level, h, fmethod = c("ets", "arima")) 
{
  data = object$y$y
  p = nrow(data)
  n = ncol(data)
  ncomp = dim(object$basis)[2] - 1
  if((n - ncomp - h + 1) <= 0)
  {
      ncomp = 1
  }
  mdata = apply(data, 1, mean)
  mdata2 = array(rep(as.matrix(mdata), B * h), dim = c(p, B, h))
  sdata = scale(t(data), scale = FALSE)
  load = as.matrix(svd(sdata)$v[, 1:ncomp])
  sco = sdata %*% load
  olivia = matrix(, ncomp, h)
  if (fmethod == "ets") {
    for (i in 1:ncomp) {
      olivia[i, ] = forecast(ets(sco[, i]), h = h)$mean
    }
  }
  if (fmethod == "arima") {
    for (i in 1:ncomp) {
      olivia[i, ] = forecast(auto.arima(sco[, i]), h = h)$mean
    }
  }
  forerr = matrix(NA, (n - ncomp - h + 1), ncomp)
  for (i in h:(n - ncomp)) {
    k = i + (ncomp - h)
    fore = matrix(, 1, ncomp)
    if (fmethod == "ets") {
      for (j in 1:ncomp) {
        fore[, j] = forecast(ets(sco[1:k, j]), h = h)$mean[h]
      }
    }
    if (fmethod == "arima") {
      for (j in 1:ncomp) {
        fore[, j] = forecast(auto.arima(sco[1:k, j]), 
                             h = h)$mean[h]
      }
    }
    forerr[i - h + 1, ] = sco[k + h, ] - fore
  }
  resi = t(sdata) - load %*% t(sco)
  q = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:p) {
      q[i, , j] = sample(resi[i, ], size = B, replace = TRUE)
    }
  }
  ny = array(NA, dim = c(ncomp, B, h))
  for (j in 1:h) {
    for (i in 1:ncomp) {
      ny[i, , j] = sample(forerr[, i], size = B, replace = TRUE)
    }
  }
  oli = array(rep(olivia, B * h), dim = c(ncomp, B, h))
  fo = array(NA, dim = c(ncomp, B, h))
  for (j in 1:h) {
    for (i in 1:B) {
      fo[, i, j] = oli[, i, j] + ny[, i, j]
    }
  }
  pred = array(NA, dim = c(p, B, h))
  for (j in 1:h) {
    for (i in 1:B) {
      pred[, i, j] = load %*% fo[, i, j] + mdata2[, i, 
                                                  j] + q[, i, j]
    }
  }
  k1 = k2 = matrix(NA, p, h)
  for (j in 1:h) {
    for (i in 1:p) {
      k1[i, j] = quantile(pred[i, , j], (100 - level)/200, 
                          na.rm = TRUE)
      k2[i, j] = quantile(pred[i, , j], 1 - (100 - level)/200, 
                          na.rm = TRUE)
    }
  }
  return(list(bootsamp = pred, lb = k1, ub = k2))
}
