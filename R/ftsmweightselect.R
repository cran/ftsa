ftsmweightselect <-
function(data, errortype = c("L1norm", "L2norm", "Linfnorm"), h)
{
  han = function(beta)
  {
    data = data$y
    n = dim(data)[2]
    fitdata = data[,1:(n-h)]
    testdata = data[,(n-h+1):n]
    s = as.numeric(rownames(fitdata))
    fore = forecast(ftsm(fts(s, fitdata), weight = TRUE, beta = beta), h = h)$mean$y
    L1norm = apply(abs(fore - testdata), 2, mean)
    L2norm = sqrt(apply((fore - testdata)^2, 2, mean))
    Linfnorm = apply(abs(fore - testdata), 2, max)
    return(switch(errortype,
                  L1norm = mean(L1norm),
                  L2norm = mean(L2norm),
                  Linfnorm = mean(Linfnorm)))
  }
  optimize(han, c(0, 1))
}

