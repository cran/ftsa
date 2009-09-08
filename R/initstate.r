initstate = function (y, trendtype, seasontype)
{
    if (seasontype != "N") {
        m <- frequency(y)
        if (length(y) >= 3 * m)
            y.d <- decompose(y, type = switch(seasontype, A = "additive",
                             M = "multiplicative"))
        else y.d <- list(seasonal = switch(seasontype, A = y -
                         mean(y), M = y/mean(y)))
        init.seas <- rev(y.d$seasonal[2:m])
        names(init.seas) <- paste("s", 0:(m - 2), sep = "")
        if (seasontype == "A")
            y.sa <- y - y.d$seasonal
        else y.sa <- y/y.d$seasonal
    }
    else {
        m <- 1
        init.seas <- NULL
        y.sa <- y
    }
    maxn <- min(max(10, 2 * m), length(y.sa))
    fit <- lsfit(1:maxn, y.sa[1:maxn])
    l0 <- fit$coef[1]
    if (trendtype == "A")
        b0 <- fit$coef[2]
    else if (trendtype == "M") {
        b0 <- 1 + fit$coef[2]/fit$coef[1]
        if (abs(b0) > 1e+10)
            b0 <- sign(b0) * 1e+10
    }
    else b0 <- NULL
    names(l0) <- "l"
    if (!is.null(b0))
        names(b0) <- "b"
    return(c(l0, b0, init.seas))
}
