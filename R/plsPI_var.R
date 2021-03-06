plsPI_var = function (data, newdata, order=2, B, alpha, lambda, var_type = "const", level = 80) 
{
  data = data$y
  n = ncol(data)
  p = nrow(data)
  n2 = length(newdata)
  newdata2 = scale(t(data), scale = FALSE)
  mdata = apply(data, 1, mean)
  mdata1 = mdata[1:n2]
  mdata2 = mdata[(n2 + 1):p]
  mdata3 = matrix(rep(mdata2, B), length(mdata2), B)
  q = matrix(, order, (n - order - 3))
  for (k in 1:(n - order - 3)) {
    j = (order + 2) + k
    load = svd(newdata2[1:j, ])$v[, 1:order]
    load2 = svd(newdata2[1:(j + 1), ])$v[, 1:order]
    sco2 = newdata2[1:(j + 1), ] %*% load2
    sco = newdata2[1:j, ] %*% load
    colnames(sco) = c("s1","s2")
    var_pred = predict(vars::VAR(sco, p = 1, type = "const"), n.ahead = 1, 
                       ci = level/100)
    fore_var = matrix(, order, 1)
    for(i in 1:order) {
      var_fit_pred = var_pred$fcst[[i]]
      fore_var[i, ] = var_fit_pred[, 1]
    }
    q[, k] = sco2[(j + 1), ] - fore_var
  }
  oldata = newdata2[1:n, ]
  mhandata = apply(data[, 1:n], 1, mean)
  mhandata = matrix(rep(mhandata, B), p, B)
  load3 = svd(oldata)$v[, 1:order]
  load4 = as.matrix(load3)[1:n2, ]
  load5 = as.matrix(load3)[(n2 + 1):p, ]
  sco3 = oldata %*% load3
  colnames(sco3) = c("s1","s2")
  var_pred = predict(vars::VAR(sco3, p = 1, type = var_type), n.ahead = 1, 
                     ci = level/100)
  fore2_var = matrix(, order, 1)
  for (i in 1:order) {
    var_fit_pred = var_pred$fcst[[i]]
    fore2_var[i, ] = var_fit_pred[, 1]
  }
  k = matrix(, order, B)
  for (i in 1:order) {
    k[i, ] = sample(q[i, ], size = B, replace = TRUE)
  }
  fore2_var = matrix(rep(fore2_var, B), order, B)
  bbeta = fore2_var + k
  resi = oldata - sco3 %*% t(load3)
  k2 = matrix(, (p - n2), B)
  for (i in 1:(p - n2)) {
    j = n2 + i
    k2[i, ] = sample(resi[, j], size = B, replace = TRUE)
  }
  I = diag(1:order)
  betapls = matrix(, order, B)
  if (n2 == 1) {
    for (i in 1:B) {
      betapls[, i] = ginv(t(matrix(load4, nrow = 1)) %*% 
                            matrix(load4, nrow = 1) + lambda * I) %*% (t(matrix(load4, 
                                                                                nrow = 1)) %*% (newdata - mdata1) + lambda * 
                                                                         bbeta[, i])
    }
  }
  if (n2 > 1) {
    for (i in 1:B) {
      betapls[, i] = ginv(t(load4) %*% load4 + lambda * 
                            I) %*% (t(load4) %*% (newdata - mdata1) + lambda * 
                                      bbeta[, i])
    }
  }
  bootsamp = load5 %*% betapls + k2 + mdata3
  w1 = w2 = w3 = w4 = matrix(, (p - n2), 1)
  for (i in 1:(p - n2)) {
    w1[i, ] = quantile(bootsamp[i, ], alpha/2)
  }
  for (i in 1:(p - n2)) {
    w2[i, ] = quantile(bootsamp[i, ], 1 - alpha/2)
  }
  forecasts = apply(bootsamp, 1, mean)
  return(list(forecasts = forecasts, bootsamp = bootsamp, low = w1, 
              up = w2))
}