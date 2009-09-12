SeasDummy = function (x)
{
    n <- length(x)
    m <- frequency(x)
    if (m == 1)
        stop("Non-seasonal data")
    tt <- 1:n
    fmat <- matrix(NA, nrow = n, ncol = 2 * m)
    for (i in 1:m) {
        fmat[, 2 * i] <- sin(2 * pi * i * tt/m)
        fmat[, 2 * (i - 1) + 1] <- cos(2 * pi * i * tt/m)
    }
    return(fmat[, 1:(m - 1)])
}
