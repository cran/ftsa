SD.test=function (wts, s = frequency(wts))
{
    if (s == 1)
        stop("Not seasonal data")
    t0 <- start(wts)
    N <- length(wts)
    if (N <= s)
        stop("Insufficient data")
    frec <- rep(1, as.integer((s + 1)/2))
    ltrunc <- round(s * (N/100)^0.25)
    R1 <- as.matrix(SeasDummy(wts))
    lmch <- lm(wts ~ R1)
    Fhat <- Fhataux <- matrix(nrow = N, ncol = s - 1)
    for (i in 1:(s - 1)) Fhataux[, i] <- R1[, i] * lmch$residuals
    for (i in 1:N) {
        for (n in 1:(s - 1)) Fhat[i, n] <- sum(Fhataux[1:i, n])
    }
    wnw <- 1 - seq(1, ltrunc, 1)/(ltrunc + 1)
    Ne <- nrow(Fhataux)
    Omnw <- 0
    for (k in 1:ltrunc) Omnw <- Omnw + (t(Fhataux)[, (k + 1):Ne] %*%
        Fhataux[1:(Ne - k), ]) * wnw[k]
    Omfhat <- (crossprod(Fhataux) + Omnw + t(Omnw)) / Ne
    sq <- seq(1, s - 1, 2)
    frecob <- rep(0, s - 1)
    for (i in 1:length(frec)) {
        if (frec[i] == 1 && i == as.integer(s/2))
            frecob[sq[i]] <- 1
        if (frec[i] == 1 && i < as.integer(s/2))
            frecob[sq[i]] <- frecob[sq[i] + 1] <- 1
    }
    a <- length(which(frecob == 1))
    A <- matrix(0, nrow = s - 1, ncol = a)
    j <- 1
    for (i in 1:(s - 1)) if (frecob[i] == 1) {
        A[i, j] <- 1
        ifelse(frecob[i] == 1, j <- j + 1, j <- j)
    }
    stL <- (1/N^2) * sum(diag(solve(t(A) %*% Omfhat %*% A, tol = 1e-25) %*%
            t(A) %*% t(Fhat) %*% Fhat %*% A))
    return(stL)
}
