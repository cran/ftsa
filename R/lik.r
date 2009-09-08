lik = function (par, y, nstate, errortype, trendtype, seasontype, damped,
                par.noopt, lowerb, upperb, opt.crit, nmse, bounds, m, 
                pnames, pnames2)
{
    names(par) <- pnames
    names(par.noopt) <- pnames2
    alpha <- c(par["alpha"], par.noopt["alpha"])["alpha"]
    if (is.na(alpha))
        stop("alpha problem!")
    if (trendtype != "N") {
        beta <- c(par["beta"], par.noopt["beta"])["beta"]
        if (is.na(beta))
            stop("beta Problem!")
    }
    else beta <- NULL
    if (seasontype != "N") {
        gamma <- c(par["gamma"], par.noopt["gamma"])["gamma"]
        if (is.na(gamma))
            stop("gamma Problem!")
    }
    else gamma <- NULL
    if (damped) {
        phi <- c(par["phi"], par.noopt["phi"])["phi"]
        if (is.na(phi))
            stop("phi Problem!")
    }
    else phi <- NULL
    if (!check.param(alpha, beta, gamma, phi, lowerb, upperb, bounds, m))
        return(1e+12)
    np <- length(par)
    init.state <- par[(np - nstate + 1):np]
    if (seasontype != "N")
        init.state <- c(init.state, m * (seasontype == "M") -
            sum(init.state[(2 + (trendtype != "N")):nstate]))
    if (seasontype == "M") {
        seas.states <- init.state[-(1:(1 + (trendtype != "N")))]
        if (min(seas.states) < 0)
            return(1e+08)
    }
    e <- pegelsresid.C(y, frequency(y), init.state, errortype,
        trendtype, seasontype, damped, alpha, beta, gamma, phi)
    if (is.na(e$lik))
        return(1e+08)
    if (e$lik < -1e+10)
        return(-1e+10)
    if (opt.crit == "lik")
        return(e$lik)
    else if (opt.crit == "mse")
        return(e$amse[1])
    else if (opt.crit == "amse")
        return(mean(e$amse[1:nmse]))
    else if (opt.crit == "sigma")
        return(mean(e$e^2))
}
