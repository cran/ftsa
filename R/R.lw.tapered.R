R.lw.tapered <-
function(d, peri, m, p, T)
{
    lambda <- 2*pi/T
    j <- seq(p,m,p)
    K <- log(p/m*(sum(peri[j]*(lambda*(j))^(2*d))))-2*d*p/m*sum(log(lambda*(j)))
    return(K)
}
