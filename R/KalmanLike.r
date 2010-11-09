KalmanLike = function (y, mod, nit = 0, fast = TRUE)
{
    x <- .Call("KalmanLike", y, mod$Z, mod$a, mod$P, mod$T, mod$V,
        mod$h, mod$Pn, as.integer(nit), FALSE, fast = fast, PACKAGE = "stats")
    names(x) <- c("ssq", "sumlog")
    s2 <- x[1]/length(y)
    list(Lik = 0.5 * (log(x[1]/length(y)) + x[2]/length(y)), s2 = s2)
}
