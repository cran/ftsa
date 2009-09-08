initparam = function (alpha, beta, gamma, phi, trendtype, seasontype, damped,
    lower, upper)
{
    par <- numeric(0)
    if (is.null(alpha)) {
        if (is.null(beta) & is.null(gamma))
            alpha <- lower[1] + 0.5 * (upper[1] - lower[1])
        else if (is.null(gamma))
            alpha <- beta + 0.001
        else if (is.null(beta))
            alpha <- 0.999 - gamma
        else alpha <- 0.5 * (beta - gamma + 1)
        if (alpha < lower[1] | alpha > upper[1])
            stop("Inconsistent parameter limits")
        par <- alpha
        names(par) <- "alpha"
    }
    if (is.null(beta)) {
        if (trendtype != "N") {
            beta <- lower[2] + 0.1 * (upper[2] - lower[2])
            if (beta > alpha)
                beta <- min(alpha - 0.001, 0.001)
            if (beta < lower[2] | beta > upper[2])
                stop("Can't find consistent starting parameters")
            par <- c(par, beta)
            names(par)[length(par)] <- "beta"
        }
    }
    if (is.null(gamma)) {
        if (seasontype != "N") {
            gamma <- lower[3] + 0.01 * (upper[3] - lower[3])
            if (gamma > 1 - alpha)
                gamma <- min(0.999 - alpha, 0.001)
            if (gamma < lower[3] | gamma > upper[3])
                stop("Can't find consistent starting parameters")
            par <- c(par, gamma)
            names(par)[length(par)] <- "gamma"
        }
    }
    if (is.null(phi)) {
        if (damped) {
            phi <- lower[4] + 0.99 * (upper[4] - lower[4])
            par <- c(par, phi)
            names(par)[length(par)] <- "phi"
        }
    }
    return(par)
}
