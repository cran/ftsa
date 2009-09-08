pegelsresid.C = function (y, m, init.state, errortype, trendtype, seasontype,
    damped, alpha, beta, gamma, phi)
{
    n <- length(y)
    p <- length(init.state)
    x <- numeric(p * (n + 1))
    x[1:p] <- init.state
    e <- numeric(n)
    lik <- 0
    if (!damped)
        phi <- 1
    if (trendtype == "N")
        beta <- 0
    if (seasontype == "N")
        gamma <- 0
    Cout <- .C("etscalc", as.double(y), as.integer(n), as.double(x),
        as.integer(m), as.integer(switch(errortype, A = 1, M = 2)),
        as.integer(switch(trendtype, N = 0, A = 1, M = 2)), as.integer(switch(seasontype,
            N = 0, A = 1, M = 2)), as.double(alpha), as.double(beta),
        as.double(gamma), as.double(phi), as.double(e), as.double(lik),
        as.double(numeric(10)), PACKAGE = "forecast")
    if (!is.na(Cout[[13]])) {
        if (abs(Cout[[13]] + 99999) < 1e-07)
            Cout[[13]] <- NA
    }
    return(list(lik = Cout[[13]], amse = Cout[[14]], e = Cout[[12]],
        states = matrix(Cout[[3]], nrow = n + 1, ncol = p, byrow = TRUE)))
}
