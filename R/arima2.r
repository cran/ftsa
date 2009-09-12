arima2 = function (x, model)
{
    use.drift <- is.element("drift", names(model$coef))
    use.intercept <- is.element("intercept", names(model$coef))
    if (use.drift)
        xreg <- length(model$x) + (1:length(x))
    if (model$arma[5] > 1 & sum(abs(model$arma[c(3, 4, 7)])) >
        0) {
        if (use.drift)
            refit <- Arima(x, order = model$arma[c(1, 6, 2)],
                seasonal = list(order = model$arma[c(3, 7, 4)],
                  period = model$arma[5]), fixed = model$coef,
                include.mean = use.intercept, xreg = xreg)
        else refit <- Arima(x, order = model$arma[c(1, 6, 2)],
            seasonal = list(order = model$arma[c(3, 7, 4)], period = model$arma[5]),
            fixed = model$coef, include.mean = use.intercept)
    }
    else if (length(model$coef) > 0) {
        if (use.drift)
            refit <- Arima(x, order = model$arma[c(1, 6, 2)],
                fixed = model$coef, xreg = xreg, include.mean = use.intercept)
        else refit <- Arima(x, order = model$arma[c(1, 6, 2)],
            fixed = model$coef, include.mean = use.intercept)
    }
    else refit <- Arima(x, order = model$arma[c(1, 6, 2)], include.mean = FALSE)
    refit$var.coef <- matrix(0, length(refit$coef), length(refit$coef))
    return(refit)
}
