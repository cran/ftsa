bin2 = function (x, ab, nbin = c(20, 20))
{
    if (missing(ab)) {
        ab <- t(array(c(nicerange(x[, 1]), nicerange(x[, 2])),
            c(2, 2)))
    }
    n <- nrow(x)
    r <- .Fortran("bin2", as.double(x), as.integer(n), as.double(ab),
        as.integer(nbin[1]), as.integer(nbin[2]), nc = integer(nbin[1] *
            nbin[2]), nskip = integer(1), PACKAGE = "ash")
    list(nc = matrix(r$nc, nbin[1], nbin[2]), ab = ab, nskip = r$nskip)
}
