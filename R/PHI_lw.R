PHI_lw <-
function(ht, T, p)
{
    k <- seq(0, T-p, p)
    lambdak <- 2*pi/T*k
    enumerator <- 0
    for(i in 1:length(k))
    {
        enumerator <- enumerator+sum(ht^2*cos((1:T)*lambdak[i]))^2
    }
    denominator <- sum(ht^2)^2
    return(enumerator/denominator)
}
