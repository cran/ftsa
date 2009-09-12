check.param = function (alpha, beta, gamma, phi, lower, upper, bounds, m)
{
    if (bounds != "admissible") {
        if (!is.null(alpha)) {
            if (alpha < lower[1] | alpha > upper[1])
                return(0)
        }
        if (!is.null(beta)) {
            if (beta < lower[2] | beta > alpha | beta > upper[2])
                return(0)
        }
        if (!is.null(phi)) {
            if (phi < lower[4] | phi > upper[4])
                return(0)
        }
        if (!is.null(gamma)) {
            if (gamma < lower[3] | gamma > 1 - alpha | gamma >
                upper[3])
                return(0)
        }
    }
    if (bounds != "usual") {
        if (!admissible(alpha, beta, gamma, phi, m))
            return(0)
    }
    return(1)
}
