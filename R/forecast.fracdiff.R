# Forecast the output of fracdiff() or arfima()

forecast.fracdiff <- function(object, h=10, level=c(80,95), fan=FALSE, ...) 
{
    # Extract data
    if (is.element("x", names(object))) 
        x <- object$x
    else 
        x <- object$x <- eval.parent(parse(text=as.character(object$call)[2]))
    n <- length(x)
    meanx <- mean(x)
    x <- x - meanx
    
    # Construct ARMA part of model and forecast with it
    y <- diffseries(x, d=object$d)
    fit <- arima(y, order=c(length(object$ar),0,length(object$ma)), include.mean=FALSE, fixed=c(object$ar,-object$ma))
    fcast.y <- forecast(fit, h=h, level=level)

    # Binomial coefficient for expansion of d
    bin.c <- (-1)^(0:(n+h)) * choose(object$d,(0:(n+h)))

    # Cumulative forecasts of y and forecast of y
    b <- numeric(n)
    fcast.x <- LHS <- numeric(h)
    RHS <- cumsum(fcast.y$mean)
    bs <- cumsum(bin.c[1:h])
    b <- bin.c[(1:n)+1]
    fcast.x[1] <- RHS[1] <- fcast.y$mean[1] - sum(b*rev(x))
    for (k in 2:h)
    {
        b <- b + bin.c[(1:n)+k]
        RHS[k] <- RHS[k] - sum(b*rev(x))
        LHS[k] <- sum(rev(fcast.x[1:(k-1)]) * bs[2:k])
        fcast.x[k] <- RHS[k] - LHS[k]
    }
    
    # Extract stuff from ARMA model
    p <- fit$arma[1]
    q <- fit$arma[2]
    phi <- theta <- numeric(h)
    if(p > 0)
        phi[1:p] <- fit$coef[1:p]
    if(q > 0)
        theta[1:q] <- fit$coef[p+(1:q)]

    # Calculate psi weights
    new.phi <- psi <- numeric(h)
    psi[1] <- new.phi[1] <- 1
    new.phi[2:h] <- -bin.c[2:h]
    for (i in 2:h) 
    {
        if(p>0)
            new.phi[i] <- sum(phi[1:(i-1)] * bin.c[(i-1):1]) - bin.c[i]
        psi[i] <- sum(new.phi[2:i] * rev(psi[1:(i-1)])) + theta[i-1]
    }
    
    # Compute forecast variances
    fse <- sqrt(cumsum(psi^2) * fit$sigma2)
    
    # Compute prediction intervals
    if (fan) 
        level <- seq(51, 99, by = 3)
    else 
    {
        if (min(level) > 0 & max(level) < 1) 
            level <- 100 * level
        else if (min(level) < 0 | max(level) > 99.99) 
            stop("Confidence limit out of range")
    }
    nint <- length(level)
    upper <- lower <- matrix(NA, ncol = nint, nrow=h)
    for (i in 1:nint) 
    {
        qq <- qnorm(0.5 * (1 + level[i]/100))
        lower[, i] <- fcast.x - qq * fse
        upper[, i] <- fcast.x + qq * fse
    }
    colnames(lower) = colnames(upper) = paste(level, "%", sep = "")

    res <- residuals(fit)
    data.tsp <- tsp(x)
    if(is.null(data.tsp))
        data.tsp <- c(1,n,1)
    mean.fcast <- ts(fcast.x+meanx, f=data.tsp[3], s=data.tsp[2] + 1/data.tsp[3])
    lower <- ts(lower+meanx, f=data.tsp[3], s=data.tsp[2] + 1/data.tsp[3])
    upper <- ts(upper+meanx, f=data.tsp[3], s=data.tsp[2] + 1/data.tsp[3])
    method <- paste("ARFIMA(",p,",",round(object$d,2),",",q,")")
    return(structure(list(x=x+meanx, mean=mean.fcast, upper=upper, lower=lower, 
        level=level, method=method, xname=deparse(substitute(x)), model=object, 
        residuals=res, fitted=x-res), class="forecast"))
}
