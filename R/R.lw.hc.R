R.lw.hc <-
function(d, peri, m, T) 
{
    j <- (1:m) + 1/2
    lambdaj <- 2 * pi *j /T
    K <- log(mean(peri[1:m] * (lambdaj )^(2 * d))) -  2 * d*mean(log(lambdaj ))
    return(K)
}
