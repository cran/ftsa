`forecast.ftsm` <- function (object, h = 10, method = c("ets", "arima", "ar", "ets.na", 
    "rwdrift", "rw"), level = 80, jumpchoice = c("fit", "actual"), pimethod = c("parametric", "nonparametric"), 
     B = 100, usedata = nrow(object$coeff), adjust = TRUE, model = NULL, damped = NULL, stationary = FALSE, ...) 
{
    method <- match.arg(method)
    jumpchoice <- match.arg(jumpchoice)
    pimethod <- match.arg(pimethod)
    if (jumpchoice == "actual") {
        var.col <- apply(object$coeff, 2, var)
        idx <- order(var.col)[1]
        if (var.col[idx] > 1e-08) 
            stop("No mean function fitted. So jumpchoice cannot be 'actual'")
        trueval <- object$fitted$y + object$residuals$y
        n <- ncol(trueval)
        object$basis[, idx] <- trueval[, n]
        for (i in 1:ncol(object$basis)) {
            if (i != idx) 
                object$coeff[, i] <- object$coeff[, i] - object$coeff[n, 
                  i]
        }
    }
    else if (jumpchoice != "fit") 
        stop("Unknown jump choice")
    nb <- ncol(object$basis)
    l <- nrow(object$coeff)
    if (method == "ar") 
        stationary <- TRUE
    if (any(stationary) == TRUE) 
        method = "arima"
    meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = nb)
    obs <- fitted <- matrix(NA, nrow = l, ncol = nb)
    qconf <- qnorm(0.5 + level/200)
    usedata <- min(usedata, l)
    ytsp <- tsp(object$coeff)
    x <- xx <- ts(as.matrix(object$coeff[(l - usedata + 1):l, 
        ]), s = ytsp[1] + l - usedata, f = ytsp[3])
    xx[object$weights[(l - usedata + 1):l] < 0.1, ] <- NA
    fmodels <- list()
    if (method == "ets") {
        if (is.null(model)) 
            model <- rep("ZZZ", nb)
        else if (length(model) == 1) 
            model <- c("ANN", rep(model, nb))
        else if (length(model) == nb - 1) 
            model <- c("ANN", model)
        else stop("Length of model does not match number of coefficients")
        if (!is.null(damped)) {
            if (length(damped) == 1) 
                damped <- c(FALSE, rep(damped, nb))
            else if (length(damped) == nb - 1) 
                damped <- c(FALSE, damped)
            else stop("Length of damped does not match number of coefficients")
        }
        for (i in 1:nb) {
            if (!is.null(damped)) 
                fmodels[[i]] <- ets(x[, i], model = model[i], 
                  damped = damped[i], ...)
            else fmodels[[i]] <- ets(x[, i], model = model[i], 
                ...)
            pegelsfit <- forecast(fmodels[[i]], h = h, level = level)
            meanfcast[, i] <- pegelsfit$mean
            varfcast[, i] <- ((pegelsfit$upper[, 1] - pegelsfit$lower[, 
                1])/(2 * qconf[1]))^2
            fitted[, i] <- pegelsfit$fitted
        }
    }
    else if (method == "ets.na") {
        for (i in 1:nb) {
            pegelsfit <- pegelsna(xx[, i], h = h, level = level, 
                ...)
            fmodels[[i]] <- pegelsfit
            meanfcast[, i] <- pegelsfit$mean
            varfcast[, i] <- ((pegelsfit$upper[, 1] - pegelsfit$lower[, 
                1])/(2 * qconf[1]))^2
            fitted[, i] <- pegelsfit$fitted
        }
    }
    else if (method == "arima") {
        if (length(stationary) == 1) 
            stationary <- c(TRUE, rep(stationary, nb))
        else if (length(stationary) == nb - 1) 
            stationary <- c(TRUE, stationary)
        else stop("Length of stationary does not match number of coefficients")
        for (i in 1:nb) {
            if (var(xx[, i], na.rm = TRUE) < 1e-08) {
                cc <- mean(xx[, i], na.rm = TRUE)
                fmodels[[i]] <- list("Constant", cc)
                meanfcast[, i] <- rep(cc, h)
                varfcast[, i] <- rep(0, h)
                fitted[, i] <- rep(cc, length(xx[, i]))
            }
            else {
                if (stationary[i]) {
                  .xxi <- xx[, i]
                  old.warn <- options(warn = -1)
                  my.ar <- function(x, max.p = 5, aic = TRUE, 
                    ...) {
                    .xxi <- x
                    fit <- ar(.xxi, method = "mle", order = max.p, 
                      aic = aic)
                    assign(".xxi", .xxi, envir = .GlobalEnv)
                    return(fit)
                  }
                  barima <- my.ar(xx[, i], ...)
                  options(old.warn)
                  fitted[, i] <- xx[, i] - barima$resid
                }
                else {
                  barima <- auto.arima(xx[, i], ...)
                  fitted[, i] <- fitted(barima)
                }
                pred <- forecast(barima, h = h, level = level)
                fmodels[[i]] <- pred
                meanfcast[, i] <- pred$mean
                varfcast[, i] <- ((pred$upper[, 1] - pred$lower[, 
                  1])/(2 * qconf[1]))^2
            }
        }
    }
    else if (method == "rwdrift") {
        for (i in 1:nb) {
            fcast <- fmodels[[i]] <- rw.drift(x[, i], ...)
            meanfcast[, i] <- object$coeff[l, i] + fcast$drift * 
                (1:h)
            varfcast[, i] <- (1:h) * fcast$see^2 + ((1:h) * fcast$sec)^2
            fitted[, i] <- fcast$fits
        }
    }
    else if (method == "rw") {
        for (i in 1:nb) {
            fmodels[[i]] <- list("Random walk", var(diff(x[, 
                i])))
            meanfcast[, i] <- rep(object$coeff[l, i], h)
            varfcast[, i] <- (1:h) * fmodels[[i]][[2]]
            fitted[, i] <- c(NA, diff(x[, i]))
        }
    }
    else stop("Unknown method")
    ytsp <- tsp(object$fitted$time)
    error <- ts(object$coeff - fitted, s = ytsp[1], f = ytsp[3])
    ferror <- onestepfcast <- object$y
    onestepfcast$y <- object$basis %*% t(fitted)
    onestepfcast$yname <- "One step forecasts"
    ferror$y <- object$y$y - onestepfcast$y
    ferror$yname <- "One step errors"
    fmean <- fts(object$y$x, object$basis %*% t(meanfcast), s = ytsp[2] + 
        1, f = ytsp[3], xname = object$y$xname, yname = "Forecasts")
    res <- object$residuals
    res$y <- res$y^2
    vx <- rowMeans(res$y)
    modelvar <- object$basis^2 %*% t(varfcast)
    totalvar <- sweep(modelvar, 1, vx + object$mean.se^2, "+")
    if (adjust & nb > 1) {
        adj.factor <- rowMeans(ferror$y^2, na.rm = TRUE)/totalvar[, 
            1]
        totalvar <- sweep(totalvar, 1, adj.factor, "*")
    }
    else adj.factor <- 1
    if (length(qconf) > 1) 
        stop("Multiple confidence levels not yet implemented")
    if (pimethod == "parametric") {
        tmp <- qconf * sqrt(totalvar)
        flower <- fts(object$y$x, fmean$y - tmp, s = ytsp[2] + 
            1, f = ytsp[3], xname = object$y$xname, yname = "Forecast lower limit")
        fupper <- fts(object$y$x, fmean$y + tmp, s = ytsp[2] + 
            1, f = ytsp[3], xname = object$y$xname, yname = "Forecast upper limit")
        coeff <- list()
        for (i in 1:nb) {
            coeff[[i]] <- structure(list(mean = ts(meanfcast[, 
                i], s = ytsp[2] + 1, f = ytsp[3]), lower = ts(meanfcast[, 
                i] - qconf * sqrt(varfcast[, i]), s = ytsp[2] + 
                1, f = ytsp[3]), upper = ts(meanfcast[, i] + 
                qconf * sqrt(varfcast[, i]), s = ytsp[2] + 1, 
                f = ytsp[3]), level = level, x = x[, i], method = method, 
                model = fmodels[[i]]), class = "forecast")
        }
        names(coeff) <- paste("Basis", 1:nb)
        return(structure(list(mean = fmean, lower = flower, upper = fupper, 
               fitted = onestepfcast, error = ferror, coeff = coeff, 
               coeff.error = error, var = list(model = modelvar, 
               error = vx, mean = object$mean.se^2, total = totalvar, 
               coeff = varfcast, adj.factor = adj.factor), model = object), 
               class = "ftsf"))
    }
    else {
         junk = ftsmPI(object, B = B, level = level, h = h, fmethod = "ets")
         lb = fts(object$y$x, junk$lb, s = ytsp[2] + 1, f = ytsp[3], xname = object$y$xname, yname = "Forecast lower limit")
         ub = fts(object$y$x, junk$ub, s = ytsp[2] + 1, f = ytsp[3], xname = object$y$xname, yname = "Forecast upper limit")
         return(structure(list(mean = fmean, bootsamp = junk$bootsamp, 
                lower = lb, upper = ub, model = object), class = "ftsf"))
    }
}
