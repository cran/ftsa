ftsm <- function (y, order = 6, ngrid = max(500, ncol(y$y)), method = c("classical", 
    "M", "rapca"), mean = TRUE, level = FALSE, lambda = 3, weight = FALSE, 
    beta = 0.1, ...) 
{
    method <- match.arg(method)
    if (!mean & !level & order < 1) 
        stop("No model to fit")
    if (order < 0) 
        stop("Order must not be negative")
    n <- ncol(y$y)
    m <- length(y$x)
    if (weight == FALSE) {
        y.pca <- fdpca(y$x, y$y, order = order, ngrid = ngrid, 
            method = method, mean = mean, level = level, lambda = lambda, 
            ...)
    }
	mean.se = fdpca(y$x, y$y, order = order, ngrid = ngrid, 
                    method = method, mean = mean, level = level, lambda = lambda, 
                    ...)$mean.se
    if (weight == TRUE) {
        newy = scale(t(y$y), scale = FALSE)
        my = apply(y$y, 1, mean)
        x = y$x
        y2 = y$y
        n = dim(y$y)[2]
        q = matrix(, n, 1)
        for (i in 1:n) {
            q[i, ] = beta * (1 - beta)^(i - 1)
        }
        wq = diag(rev(q))
        newy2 = wq %*% newy
        load = svd(newy2)$v[,1:(order)]
        sco = newy%*%load
    }
    ytsp <- tsp(y$time)
    if (weight == TRUE) {
        coeff <- ts(cbind(rep(1, dim(y$y)[2]), sco), 
                    s = ytsp[1], f = ytsp[3])
    }
    else {
        coeff <- ts(y.pca$coeff, s = ytsp[1], f = ytsp[3])
    }
    if (weight == TRUE) {
        basis <- cbind(my, load)
        colnames(basis) = c("mean", paste("phi", 1:order, sep=""))
        fits <- fts(y$x, load %*% t(sco) + my, s = ytsp[1], 
            f = ytsp[3], xname = y$xname, yname = paste("Fitted", 
                y$yname))
    }
    else {
        basis <- y.pca$basis
        fits <- fts(y$x, basis %*% t(coeff), s = ytsp[1], f = ytsp[3], 
            xname = y$xname, yname = paste("Fitted", y$yname))
    }
    rownames(basis) <- paste(y$x)
    res <- fts(y$x, y$y - fits$y, s = ytsp[1], f = ytsp[3], xname = y$xname, 
        yname = paste("Residuals", y$yname))
    if (weight == FALSE) {
        out <- list(x1 = as.numeric(colnames(y$y)), y1 = as.numeric(rownames(y$y)),
                    y = fts(y$x, y$y, xname=y$xname, yname=y$yname), basis = basis, coeff = coeff, 
		    fitted = fits, residuals = res, varprop = y.pca$varprop, wt = ts(y.pca$weights, 
                    s = ytsp[1], f = ytsp[3]), v = ts(y.pca$v, s = ytsp[1], 
                    f = ytsp[3]), basis2 = y.pca$basis2, coeff2 = y.pca$coeff2, 
                    mean.se = mean.se, call = match.call())
    }
    else {
        out <- list(x1 = as.numeric(colnames(y$y)), y1 = as.numeric(rownames(y$y)),
             y = fts(y$x, y2, xname=y$xname, yname=y$yname), basis = basis, coeff = coeff, 
             fitted = fits, residuals = res, wt = rev(q), mean.se = mean.se,
			 call = match.call())
    }
    return(structure(out, class = c("ftsm", "fm")))
}
