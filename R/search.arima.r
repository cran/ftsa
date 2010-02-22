search.arima = function (x, d = NA, D = NA, max.p = 5, max.q = 5, max.P = 2,
    max.Q = 2, max.order = 5, stationary = FALSE, ic = c("aic",
     "aicc", "bic"), trace = FALSE, approximation = FALSE)
{
    ic <- match.arg(ic)
    m <- frequency(x)
    seasonal <- (m > 1)
    if (stationary)
        d <- D <- 0
    if (!seasonal)
        D <- max.P <- max.Q <- 0
    else if (is.na(D))
        D <- nsdiffs(x)
    if (D > 0)
        dx <- diff(x, D)
    else dx <- x
    if (is.na(d))
        d <- ndiffs(dx)
    if (m > 1) {
        if (max.P > 0)
            max.p <- min(max.p, m - 1)
        if (max.Q > 0)
            max.q <- min(max.q, m - 1)
    }
    if (approximation) {
        if (D == 0)
            fit <- try(arima(x, order = c(1, d, 0)))
        else fit <- try(arima(x, order = c(1, d, 0), seasonal = list(order = c(0,
                        D, 0), period = m)))
        if (class(fit) != "try-error")
            offset <- -2 * fit$loglik - length(x) * log(fit$sigma2)
        else {
             warning("Unable to calculate AIC offset")
             offset <- 0
        }
    }
    else offset <- 0
    oldwarn <- options()$warn
    options(warn = -1)
    on.exit(options(warn = oldwarn))
    best.ic <- 1e+20
    for (i in 0:max.p) {
        for (j in 0:max.q) {
            for (I in 0:max.P) {
                for (J in 0:max.Q) {
                  if (i + j + I + J <= max.order) {
                    for (K in 0:(d + D <= 1)) {
                      fit <- myarima(x, order = c(i, d, j), seasonal = c(I,
                        D, J), constant = (K == 1), trace = trace,
                        ic = ic, approximation = approximation,
                        offset = offset)
                      if (fit$ic < best.ic) {
                          best.ic <- fit$ic
                          bestfit <- fit
                      }
                    }
                  }
                }
            }
        }
    }
    if (exists("bestfit")) {
        bestfit$x <- x
        if (approximation) {
            constant <- length(bestfit$coef) > sum(bestfit$arma[1:4])
            newbestfit <- myarima(x, order = bestfit$arma[c(1,
                6, 2)], seasonal = bestfit$arma[c(3, 7, 4)],
                constant = constant, ic, trace = FALSE, approximation = FALSE)
            if (newbestfit$ic > 1e+19) {
                options(warn = oldwarn)
                warning("Unable to fit final model using maximum likelihood. AIC value approximated")
            }
            else bestfit <- newbestfit
        }
    }
    else stop("No ARIMA model able to be estimated")
    bestfit$series <- deparse(substitute(x))
    bestfit$ic = NULL
    if (trace)
        cat("\n\n")
    return(bestfit)
}
