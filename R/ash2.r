ash2 = function (bins, m = c(5, 5), kopt = c(2, 2))
{
    nc <- bins$nc
    if (!is.matrix(nc))
        stop("bins does not contain bin count matrix")
    ab <- bins$ab
    if (!is.matrix(ab))
        stop("ab not a matrix - should be 2 by 2")
    nbin <- dim(nc)
    r <- .Fortran("ash2", as.integer(m[1]), as.integer(m[2]),
        as.integer(nc), as.integer(nbin[1]), as.integer(nbin[2]),
        as.double(ab), as.integer(kopt), f = double(nbin[1] *
            nbin[2]), double(m[1] * m[2]), ier = double(1), PACKAGE = "ash")
    if (r$ier == 1)
        print(" estimate nonzero outside ab rectangle")
    list(z = matrix(r$f, nbin[1], nbin[2]), x = center(ab[1,
        ], nbin[1])[[1]], y = center(ab[2, ], nbin[2])[[1]],
        ab = ab, m = m, kopt = kopt, ier = r$ier)
}
