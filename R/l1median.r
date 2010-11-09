l1median = function (X, MaxStep = 200, ItTol = 10^-8)
{
    if (class(X) != "matrix") {
        if (class(X) == "data.frame")
            X = as.matrix(X)
        else X = matrix(X, ncol = 1)
    }
    ret = .C("l1median", PACKAGE = "pcaPP", as.double(X), as.integer(nrow(X)),
        as.integer(ncol(X)), med = double(ncol(X)), ret = integer(1),
        as.integer(MaxStep), as.double(ItTol))
    if (ret$ret != 0)
        return(ret$med)
    stop("iteration failed")
}
