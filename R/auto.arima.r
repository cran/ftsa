auto.arima = function (x, d = NA, D = NA, max.p = 5, max.q = 5, max.P = 2,
    max.Q = 2, max.order = 5, start.p = 2, start.q = 2, start.P = 1,
    start.Q = 1, stationary = FALSE, ic = c("aic", "aicc", "bic"),
    stepwise = TRUE, trace = FALSE, approximation = length(x) >
        100 | frequency(x) > 12)
{
    if (!stepwise)
        return(search.arima(x, d, D, max.p, max.q, max.P, max.Q,
            max.order, stationary, ic, trace, approximation))
    ic <- match.arg(ic)
    m <- frequency(x)
    if (stationary)
        d <- D <- 0
    if (m == 1)
        D <- max.P <- max.Q <- 0
    else if (is.na(D))
        D <- forecast:::nsdiffs(x)
    if (D > 0)
        dx <- diff(x, D)
    else dx <- x
    if (is.na(d))
        d <- forecast:::ndiffs(dx)
    if (m > 1) {
        if (max.P > 0)
            max.p <- min(max.p, m - 1)
        if (max.Q > 0)
            max.q <- min(max.q, m - 1)
    }
    p <- start.p <- min(start.p, max.p)
    q <- start.q <- min(start.q, max.q)
    P <- start.P <- min(start.P, max.P)
    Q <- start.Q <- min(start.Q, max.Q)
    constant <- (d + D <= 1)
    results <- matrix(NA, nrow = 100, ncol = 8)
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
    bestfit <- myarima(x, order = c(p, d, q), seasonal = c(P,
        D, Q), constant = constant, ic, trace, approximation,
        offset = offset)
    results[1, ] <- c(p, d, q, P, D, Q, constant, bestfit$ic)
    fit <- myarima(x, order = c(0, d, 0), seasonal = c(0, D,
        0), constant = constant, ic, trace, approximation, offset = offset)
    results[2, ] <- c(0, d, 0, 0, D, 0, constant, fit$ic)
    if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- q <- P <- Q <- 0
    }
    if (max.p > 0 | max.P > 0) {
        fit <- myarima(x, order = c(max.p > 0, d, 0), seasonal = c((m >
            1) & (max.P > 0), D, 0), constant = constant, ic,
            trace, approximation, offset = offset)
        results[3, ] <- c(1, d, 0, m > 1, D, 0, constant, fit$ic)
        if (fit$ic < bestfit$ic) {
            bestfit <- fit
            p <- (max.p > 0)
            P <- (m > 1) & (max.P > 0)
            q <- Q <- 0
        }
    }
    if (max.q > 0 | max.Q > 0) {
        fit <- myarima(x, order = c(0, d, max.q > 0), seasonal = c(0,
            D, (m > 1) & (max.Q > 0)), constant = constant, ic,
            trace, approximation, offset = offset)
        results[4, ] <- c(0, d, 1, 0, D, m > 1, constant, fit$ic)
        if (fit$ic < bestfit$ic) {
            bestfit <- fit
            p <- P <- 0
            Q <- (m > 1) & (max.Q > 0)
            q <- (max.q > 0)
        }
    }
    startk <- 0
    k <- 4
    while (startk < k & k < 94) {
        startk <- k
        if (P > 0 & newmodel(p, d, q, P - 1, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                1, D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q, P - 1, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                P <- (P - 1)
            }
        }
        if (P < max.P & newmodel(p, d, q, P + 1, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                1, D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q, P + 1, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                P <- (P + 1)
            }
        }
        if (Q > 0 & newmodel(p, d, q, P, D, Q - 1, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                D, Q - 1), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q, P, D, Q - 1, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                Q <- (Q - 1)
            }
        }
        if (Q < max.Q & newmodel(p, d, q, P, D, Q + 1, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                D, Q + 1), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q, P, D, Q + 1, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                Q <- (Q + 1)
            }
        }
        if (Q > 0 & P > 0 & newmodel(p, d, q, P - 1, D, Q - 1,
            constant, results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P -
                1, D, Q - 1), constant = constant, ic, trace,
                approximation, offset = offset)
            results[k, ] <- c(p, d, q, P - 1, D, Q - 1, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                Q <- (Q - 1)
                P <- (P - 1)
            }
        }
        if (Q < max.Q & P < max.P & newmodel(p, d, q, P + 1,
            D, Q + 1, constant, results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P +
                1, D, Q + 1), constant = constant, ic, trace,
                approximation, offset = offset)
            results[k, ] <- c(p, d, q, P + 1, D, Q + 1, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                Q <- (Q + 1)
                P <- (P + 1)
            }
        }
        if (p > 0 & newmodel(p - 1, d, q, P, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p - 1, d, q), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p - 1, d, q, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                p <- (p - 1)
            }
        }
        if (p < max.p & newmodel(p + 1, d, q, P, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p + 1, d, q), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p + 1, d, q, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                p <- (p + 1)
            }
        }
        if (q > 0 & newmodel(p, d, q - 1, P, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q - 1), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q - 1, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                q <- (q - 1)
            }
        }
        if (q < max.q & newmodel(p, d, q + 1, P, D, Q, constant,
            results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q + 1), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q + 1, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                q <- (q + 1)
            }
        }
        if (q > 0 & p > 0 & newmodel(p - 1, d, q - 1, P, D, Q,
            constant, results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p - 1, d, q - 1), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p - 1, d, q - 1, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                q <- (q - 1)
                p <- (p - 1)
            }
        }
        if (q < max.q & p < max.p & newmodel(p + 1, d, q + 1,
            P, D, Q, constant, results[1:k, ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p + 1, d, q + 1), seasonal = c(P,
                D, Q), constant = constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p + 1, d, q + 1, P, D, Q, constant,
                fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                q <- (q + 1)
                p <- (p + 1)
            }
        }
        if (newmodel(p, d, q, P, D, Q, !constant, results[1:k,
            ])) {
            k <- k + 1
            fit <- myarima(x, order = c(p, d, q), seasonal = c(P,
                D, Q), constant = !constant, ic, trace, approximation,
                offset = offset)
            results[k, ] <- c(p, d, q, P, D, Q, !constant, fit$ic)
            if (fit$ic < bestfit$ic) {
                bestfit <- fit
                constant <- !constant
            }
        }
    }
    if (approximation) {
        newbestfit <- myarima(x, order = bestfit$arma[c(1, 6,
            2)], seasonal = bestfit$arma[c(3, 7, 4)], constant = constant,
            ic, trace = FALSE, approximation = FALSE)
        if (newbestfit$ic > 1e+19) {
            options(warn = oldwarn)
            warning("Unable to fit final model using maximum likelihood. AIC value approximated")
        }
        else bestfit <- newbestfit
    }
    bestfit$x <- x
    bestfit$series <- deparse(substitute(x))
    bestfit$ic = NULL
    if (trace)
        cat("\n\n Best model:", arima.string(bestfit), "\n\n")
    return(bestfit)
}
