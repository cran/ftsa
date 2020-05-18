ER_GR <-
function(data)
{
    n = nrow(data)
    p = ncol(data)
    m = min(n, p)

    covar = (t(data) %*% data)/(n * p)
    eigen_values = eigen(covar)$values

    V_0 = sum(eigen_values[1:m])
    kmax_star = length(which(eigen_values >= V_0/m))
    kmax = min(kmax_star, 0.1 * m)

    ER = vector("numeric", kmax - 1)
    for(k in 1:(kmax - 1))
    {
        ER[k] = eigen_values[k]/eigen_values[k+1]
    }
    k_ER = which.max(ER)

    GR = vector("numeric", kmax - 1)
    for(k in 1:(kmax - 1))
    {
        V_1 = sum(eigen_values[(k+1):m])
        u_star_1 = eigen_values[k]/V_1

        V_2 = sum(eigen_values[(k+2):m])
        u_star_2 = eigen_values[(k+1)]/V_2

        GR[k] = log(1 + u_star_1)/log(1 + u_star_2)
    }
    k_GR = which.max(GR)

    return(list(k_ER = k_ER, k_GR = k_GR))
}
