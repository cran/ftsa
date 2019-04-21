FPE.trace<-function(res, p = 2)
{
    d = length(res[1,])
    n = length(res[,1])
    if(d == 1)
    {
        out = (p*d+n)/(n-p*d)*var(res)
    }
    else
    {
        out = (p*d+n)/(n-p*d)*sum(diag(cov(res)))
    }
    return(out)
}

