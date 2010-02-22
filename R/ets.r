ets = function (y, model = "ZZZ", damped = NULL, alpha = NULL, beta = NULL, gamma = NULL, phi = NULL, 
       additive.only = FALSE, lower = c(rep(0.01, 3), 0.8), upper = c(rep(0.99, 3), 0.98), 
        opt.crit = c("lik", "amse", "mse", "sigma"), nmse = 3, bounds = c("both",
        "usual", "admissible"), ic = c("aic", "aicc", "bic"), restrict = TRUE)
{
    opt.crit <- match.arg(opt.crit)
    bounds <- match.arg(bounds)
    ic <- match.arg(ic)
    if (max(y) > 1e+06)
        warning("Very large numbers which may cause numerical problems. Try scaling the data first")
    if (class(y) == "data.frame" | class(y) == "list" | class(y) ==
        "matrix")
        stop("y should be a vector")
    y <- as.ts(y)
    if (nmse < 1 | nmse > 10)
        stop("nmse out of range")
    m <- frequency(y)
    if (m > 24)
        stop("Frequency too high")
    if (sum((upper - lower) > 0) < 4)
        stop("Lower limits must be less than upper limits")
    if (class(model) == "ets") {
        alpha = model$par["alpha"]
        beta = model$par["beta"]
        if (is.na(beta))
            beta <- NULL
        gamma = model$par["gamma"]
        if (is.na(gamma))
            gamma <- NULL
        phi = model$par["phi"]
        if (is.na(phi))
            phi <- NULL
        damped = (model$components[4] == "TRUE")
        model = paste(model$components[1], model$components[2],
            model$components[3], sep = "")
    }
    errortype <- substr(model, 1, 1)
    trendtype <- substr(model, 2, 2)
    seasontype <- substr(model, 3, 3)
    if (!is.element(errortype, c("M", "A", "Z")))
        stop("Invalid error type")
    if (!is.element(trendtype, c("N", "A", "M", "Z")))
        stop("Invalid trend type")
    if (!is.element(seasontype, c("N", "A", "M", "Z")))
        stop("Invalid season type")
    if (m == 1) {
        if (seasontype == "A" | seasontype == "M")
            stop("Nonseasonal data")
        else substr(model, 3, 3) <- seasontype <- "N"
    }
    if (restrict) {
        if ((errortype == "A" & (trendtype == "M" | seasontype ==
            "M")) | (errortype == "M" & trendtype == "M" & seasontype ==
            "A") | (additive.only & (errortype == "M" | trendtype ==
            "M" | seasontype == "M")))
            stop("Forbidden model combination")
    }
    data.positive <- (min(y) > 0)
    if (!data.positive & errortype == "M")
        stop("Inappropriate model for data with negative or zero values")
    if (!is.null(damped)) {
        if (damped & trendtype == "N")
            stop("Forbidden model combination")
    }
    if (errortype == "Z")
        errortype <- c("A", "M")
    if (trendtype == "Z")
        trendtype <- c("N", "A", "M")
    if (seasontype == "Z")
        seasontype <- c("N", "A", "M")
    if (is.null(damped))
        damped <- c(TRUE, FALSE)
    best.ic <- Inf
    for (i in 1:length(errortype)) {
        for (j in 1:length(trendtype)) {
            for (k in 1:length(seasontype)) {
                for (l in 1:length(damped)) {
                  if (trendtype[j] == "N" & damped[l])
                    next
                  if (restrict) {
                    if (errortype[i] == "A" & (trendtype[j] ==
                      "M" | seasontype[k] == "M"))
                      next
                    if (errortype[i] == "M" & trendtype[j] ==
                      "M" & seasontype[k] == "A")
                      next
                    if (additive.only & (errortype[i] == "M" |
                      trendtype[j] == "M" | seasontype[k] ==
                      "M"))
                      next
                  }
                  if (!data.positive & errortype[i] == "M")
                    next
                  fit <- etsmodel(y, errortype[i], trendtype[j],
                    seasontype[k], damped[l], alpha, beta, gamma,
                    phi, lower = lower, upper = upper, opt.crit = opt.crit,
                    nmse = nmse, bounds = bounds)
                  fit.ic <- switch(ic, aic = fit$aic, bic = fit$bic,
                    aicc = fit$aicc)
                  if (!is.na(fit.ic)) {
                    if (fit.ic < best.ic) {
                      model <- fit
                      best.ic <- fit.ic
                      best.e <- errortype[i]
                      best.t <- trendtype[j]
                      best.s <- seasontype[k]
                      best.d <- damped[l]
                    }
                  }
                }
            }
        }
    }
    if (best.ic == Inf)
        stop("No model able to be fitted")
    model$m <- m
    model$method <- paste("ETS(", best.e, ",", best.t, ifelse(best.d,
        "d", ""), ",", best.s, ")", sep = "")
    model$components <- c(best.e, best.t, best.s, best.d)
    model$call <- match.call()
    model$initstate <- model$states[1, ]
    model$sigma2 <- var(model$residuals, na.rm = TRUE)
    model$x <- as.ts(y)
    return(structure(model, class = "ets"))
}
