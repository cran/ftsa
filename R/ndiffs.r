ndiffs = function (x, alpha = 0.05)
{
    require(tseries)
    x <- c(na.omit(c(x)))
    d <- 0
    oldwarn <- options(warn = -1)
    p.v <- kpss.test(x)$p.value
    if (is.na(p.v)) {
        options(warn = oldwarn$warn)
        return(d)
    }
    while (p.v < alpha & d < 2) {
        x <- diff(x)
        d <- d + 1
        p.v <- kpss.test(x)$p.value
        if (is.na(p.v))
            return(d - 1)
    }
    options(warn = oldwarn$warn)
    return(d)
}
