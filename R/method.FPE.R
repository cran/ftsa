method.FPE <- function(object, D = 21, var_type = "const", Pmax)
{
    # compute the total variance
    
    eigen_values = eigen(var(t(object$y)))$values
    vartot = sum(eigen_values)
    
    # the matrix below will conain the different (FPE + vartot - var.explain) values
    # for p in 0:Pmax and d in 1:D
    
    pca.FTS = ftsm(object, order = D)
    pca.FTS_scores = pca.FTS$coeff[,2:(D+1)]
    
    values = matrix(NA, D, (Pmax + 1))
    for(d in 1:D)
    {
        scores = pca.FTS_scores[,1:d]
        var.explain = sum(eigen_values[1:d])
        for(p in 0:Pmax)
        {
            if(d == 1)
            {
                res = as.matrix(arima(scores, order = c(p, 0, 0), method = "ML")$residuals)
            }
            else
            {
                if(p == 0)
                {
                    mean = t(matrix(rep(colMeans(scores), nrow(scores)), d))
                    res = scores - mean
                }
                else
                {
                    colnames(scores) <- as.character(seq(1:d))
                    res = resid(vars::VAR(scores, p = p, type = var_type))
                }
            }
            values[d, p+1] = FPE.trace(res = res, p = p) + vartot - var.explain
        }
    }
    rownames(values) = 1:D
    colnames(values) = 0:Pmax
    
    # compute the estimates hat.p and hat.d for optimal order p and dimension d
    
    hat.p = (which.min(values) - 1) %/% D
    hat.d = which.min(values) %% D
    if(hat.d == 0) hat.d = hat.d + D
    out = c(hat.p, hat.d)
    return(out)
}
