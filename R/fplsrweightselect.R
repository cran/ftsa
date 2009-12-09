fplsrweightselect <-
function(data, errortype = c("L1norm", "L2norm", "Linfnorm"), h)
{
  han = function(beta)
  {
    data = data$y
    n = dim(data)[2]
    fitdata = data[,1:(n-h)]
    testdata = data[,(n-h+1):n]
    s = as.numeric(rownames(fitdata))
    r = dim(fitdata)[2]
    fore = matrix(, dim(fitdata)[1], h)
    fitdata = fts(s, fitdata)
    for(i in 1:h){
        fore[,i] = fplsr(fitdata, weight = TRUE, beta = beta)$Ypred$y
        fitdata = cbind(fitdata$y, fore[,i])
        colnames(fitdata) = as.numeric(colnames(fitdata)[1]):(as.numeric(colnames(fitdata)[r]) + i)
        fitdata = fts(s, fitdata)
    }             
    L1norm = apply(abs(fore - testdata), 2, mean)
    L2norm = sqrt(apply((fore - testdata)^2, 2, mean))
    Linfnorm = apply(abs(fore - testdata), 2, max)
    return(switch(errortype,
                  L1norm = mean(L1norm),
                  L2norm = mean(L2norm),
                  Linfnorm = mean(Linfnorm)))
  }
  optimize(han, c(0,1))
}

