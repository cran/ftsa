admissible = function (alpha, beta, gamma, phi, m)
{
    if (is.null(phi))
        phi <- 1
    if (phi < 0 | phi > 1 + 1e-08)
        return(0)
    if (is.null(gamma)) {
        if (alpha < 1 - 1/phi | alpha > 1 + 1/phi)
            return(0)
        if (!is.null(beta)) {
            if (beta < alpha * (phi - 1) | beta > (1 + phi) *
                (2 - alpha))
                return(0)
        }
    }
    else {
        if (is.null(beta))
            beta <- 0
        if (gamma < max(1 - 1/phi - alpha, 0) | gamma > 1 + 1/phi -
            alpha)
            return(0)
        if (alpha < 1 - 1/phi - gamma * (1 - m + phi + phi *
            m)/(2 * phi * m))
            return(0)
        if (beta < -(1 - phi) * (gamma/m + alpha))
            return(0)
        P <- c(phi * (1 - alpha - gamma), alpha + beta - alpha *
            phi + gamma - 1, rep(alpha + beta - alpha * phi,
            m - 2), (alpha + beta - phi), 1)
        roots <- polyroot(P)
        if (max(abs(roots)) > 1 + 1e-10)
            return(0)
    }
    return(1)
}
