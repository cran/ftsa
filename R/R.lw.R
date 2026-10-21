R.lw <-
function(d, peri, m, T, l = 1)
{
    lambda<-2*pi/T
    K <- log(1/(m-l-1)*(sum(peri[l:m]*(lambda*(l:m))^(2*d))))-2*d/(m-l-1)*sum(log(lambda*(l:m)))
    return(K)
}
