ftsmPI <- function (object, B, level, h, fmethod = c("ets", "arima")) 
{
    data = object$y$y
    p = dim(data)[1]
    n = dim(data)[2]
    ncomp = dim(object$basis)[2] - 1
    mdata = apply(data, 1, mean)
    mdata2 = array(rep(as.matrix(mdata), B*h), dim = c(p, B, h))
    sdata = scale(t(data), scale = FALSE)
    load = svd(sdata)$v[, 1:ncomp]
    sco = sdata %*% load
    oli = matrix(, ncomp, h)
    for (i in 1:ncomp) {
         oli[i, ] = forecast(ets(sco[, i]), h = h)$mean
    }     
    forerr = matrix(, (n - ncomp), ncomp)
    for (i in 1:(n - ncomp)) {
         k = i + (ncomp - 1)
         fore = matrix(, 1, ncomp)
         if (fmethod == "ets") {
             for (j in 1:ncomp) {
                  fore[, j] = forecast(ets(sco[1:k, j]), h = 1)$mean
             }
         }
         else {
              for (j in 1:ncomp) {
                   fore[, j] = forecast(auto.arima(sco[1:k, j]), 
                                        h = 1)$mean
              }
         }
         forerr[i, ] = sco[k + 1, ] - fore
    }
    resi = t(sdata) - load %*% t(sco)
    q = array(, dim = c(p, B, h))
    for (j in 1:h){
         for (i in 1:p) {
             q[i,,j] = sample(resi[i, ], size = B, replace = TRUE)
         }
    }
    ny = array(, dim = c(ncomp, B, h))
    for (j in 1:h){
         for (i in 1:ncomp) {
              ny[i,,j] = sample(forerr[, i], size = B, replace = TRUE)
         }
    }
    oli = array(rep(oli, B*h), dim = c(ncomp, B, h))
    fo = array(, dim = c(ncomp, B, h))
    for(j in 1:h){
        for (i in 1:B) {
             fo[, i, j] = oli[, i, j] + ny[, i, j]
        }
    }  
    pred = array(, dim = c(p, B, h))
    for(j in 1:h){
        for(i in 1:B){
            pred[,i,j] = load %*% fo[,i,j] + mdata2[,i,j] + q[,i,j]
        }
    }    
    k = k1 = matrix(, p, h)
    for(j in 1:h){
        for (i in 1:p) {
             k[i, j] = quantile(pred[i,,j], 1 - (100 - level)/200, na.rm = TRUE)
             k1[i, j] = quantile(pred[i,,j], (100 - level)/200, na.rm = TRUE)
        }
    }
    return(list(bootsamp = pred, lb = k1, ub = k))
}
  
