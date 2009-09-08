myarima = function (x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant = TRUE,
    ic = "aic", trace = FALSE, approximation = FALSE, offset = 0)
{
    n <- length(x)
    m <- frequency(x)
    use.season <- (sum(seasonal) > 0) & m > 0
    diffs <- order[2] + seasonal[2]
    if (approximation)
        method <- "CSS"
    else method <- "CSS-ML"
    if (diffs == 1 & constant) {
        xreg <- 1:length(x)
        if (use.season)
            fit <- try(stats:::arima(x = x, order = order, seasonal = list(order = seasonal,
                period = m), xreg = xreg, method = method), silent = TRUE)
        else fit <- try(stats:::arima(x = x, order = order, xreg = xreg,
            method = method), silent = TRUE)
    }
    else {
        if (use.season)
            fit <- try(stats:::arima(x = x, order = order, seasonal = list(order = seasonal,
                period = m), include.mean = constant, method = method),
                silent = TRUE)
        else fit <- try(stats:::arima(x = x, order = order, include.mean = constant,
            method = method), silent = TRUE)
    }
    if (class(fit) != "try-error") {
        if (diffs == 1 & constant) {
            fitnames <- names(fit$coef)
            fitnames[length(fitnames)] <- "drift"
            names(fit$coef) <- fitnames
            fit$call$xreg <- xreg
        }
        npar <- length(fit$coef) + 1
        if (approximation)
            fit$aic <- offset + n * log(fit$sigma2) + 2 * npar
        if (!is.na(fit$aic)) {
            fit$bic <- fit$aic + npar * (log(n) - 2)
            fit$aicc <- fit$aic + 2 * npar * (n/(n - npar - 1) -
                1)
            fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic,
                aicc = fit$aicc)
        }
        else fit$aic <- fit$bic <- fit$aicc <- fit$ic <- 1e+20
        minroot <- 2
        if (order[1] + seasonal[1] > 0) {
            testvec <- fit$model$phi
            last.nonzero <- max((1:length(testvec))[abs(testvec) > 1e-08])
            testvec <- testvec[1:last.nonzero]
            if (last.nonzero > 48)
                warning("Unable to check for unit roots")
            else minroot <- min(minroot, abs(polyroot(c(1, -testvec))))
        }
        if (order[3] + seasonal[3] > 0) {
            testvec <- fit$model$theta
            last.nonzero <- max((1:length(testvec))[abs(testvec) >
                1e-08])
            testvec <- testvec[1:last.nonzero]
            if (last.nonzero > 48)
                warning("Unable to check for unit roots")
            else minroot <- min(minroot, abs(polyroot(c(1, testvec))))
        }
        if (minroot < 1 + 0.001)
            fit$ic <- 1e+20
        if (trace)
            cat("\n", arima.string(fit), ":", fit$ic)
        return(fit)
    }
    else {
        if (trace) {
            cat("\n ARIMA(", order[1], ",", order[2], ",", order[3],
                ")", sep = "")
            if (use.season)
                cat("(", seasonal[1], ",", seasonal[2], ",",
                  seasonal[3], ")[", m, "]", sep = "")
            if (constant & (order[2] + seasonal[2] == 0))
                cat(" with non-zero mean")
            else if (constant & (order[2] + seasonal[2] == 1))
                cat(" with drift        ")
            else if (!constant & (order[2] + seasonal[2] == 0))
                cat(" with zero mean    ")
            else cat("         ")
            cat(" :", 1e+20, "*")
        }
        return(list(ic = 1e+20))
    }
}
