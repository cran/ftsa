fdpca = function (x, y, order = 2, ngrid = 500, method = "M", mean = mean,
                  level = level, lambda = 2.3262, iter = 1, ...)
{
    n <- ncol(y)
    m <- length(x)
    if (lambda < 1)
        stop("Lambda too small")
    if (iter < 1)
        stop("Need at least one iteration")
    if (order < 0)
        stop("Order must be at least 0")
    if (ngrid < n)
        stop("Grid should be larger than number of observations per time period.")
    if (m != nrow(y))
        stop("x and y of incompatible dimension")
    if (order > n/2 & method != "classical") {
        warning("Not enough data for robust PCA")
        method <- "classical"
    }
    yy <- matrix(NA, nrow = ngrid, ncol = n)
    xx <- seq(min(x), max(x), l = ngrid)
    delta <- xx[2] - xx[1]
    for (i in 1:n) {
        miss <- is.na(y[, i])
        yy[, i] <- spline(x[!miss], y[!miss, i], n = ngrid)$y
    }
    if (mean) {
        if (method == "M" | method == "rapca")
            ax <- L1median2(t(yy), method = "hoss")
        else ax <- rowMeans(yy, na.rm = TRUE)
        yy <- sweep(yy, 1, ax)
        axse <- approx(xx, sqrt(apply(yy, 1, var)/n), xout = x)$y
    }
    else axse <- NULL
    if (level) {
        bx <- colMeans(yy, na.rm = TRUE)
        yy <- sweep(yy, 2, bx)
    }
    if (level) {
        coeff <- as.matrix(bx)
        basis <- as.matrix(rep(1, m))
        colnames(coeff) <- colnames(basis) <- "level"
    }
    else coeff <- basis <- NULL
    if (mean) {
        coeff <- cbind(rep(1, n), coeff)
        basis <- cbind(approx(xx, ax, xout = x)$y, basis)
        colnames(coeff)[1] <- colnames(basis)[1] <- "mean"
    }
    if (order == 0)
        return(list(basis = basis, coeff = coeff, weights = rep(1,
            n), v = rep(1, n), mean.se = axse))
    if (method == "M" | method == "rapca") {
        robusteig <- rapca(yy, order = order, mean = FALSE, ...)
        B <- robusteig$coeff
        Phi <- robusteig$basis
        yyhat <- Phi %*% t(B)
        v <- colSums((yy - yyhat)^2) * delta
    }
    else {
        v <- rep(1, n)
        iter <- 1
    }
    if (method == "rapca") {
        varprop <- rep(NA, order)
        w <- rep(1, order)
    }
    else {
        for (i in 1:iter) {
            medv <- median(v)
            w <- ifelse(v < medv + lambda * sqrt(medv), 1, 0)
            V <- repmat(w/mean(w), 1, ngrid) * t(yy)
            s <- La.svd(V)
            Phi <- as.matrix(t(s$vt)[, 1:order])
            B <- t(yy) %*% Phi
            v <- colSums((yy - Phi %*% t(B))^2) * delta
        }
        varprop <- s$d^2
        varprop <- varprop/sum(s$d^2)
    }
    colnames(B) <- paste("beta", 1:order, sep = "")
    coeff <- cbind(coeff, B * delta)
    m <- ncol(basis)
    for (i in 1:order) {
        basis <- cbind(basis, approx(xx, Phi[, i], xout = x)$y/delta)
        if (sum(basis[, i + m]) < 0) {
            basis[, i + m] <- -basis[, i + m]
            coeff[, i + m] <- -coeff[, i + m]
        }
    }
    colnames(basis)[m + (1:order)] <- paste("phi", 1:order, sep = "")
    varprop <- varprop[1:order]
    yy <- yy - Phi %*% t(B)
    s <- try(La.svd(t(yy)), silent = TRUE)
    if (class(s) == "try-error") {
        s <- svd(t(yy), LINPACK = TRUE)
        s$vt <- t(s$v)
    }
    Phi2 <- as.matrix(t(s$vt)[, s$d > 1e-06])
    m <- ncol(Phi2)
    basis2 <- coeff2 <- NULL
    if (m > 0) {
        B2 <- t(yy) %*% Phi2
        colnames(B2) <- paste("beta", order + (1:ncol(B2)), sep = "")
        coeff2 <- B2 * delta
        for (i in 1:m) {
            basis2 <- cbind(basis2, approx(xx, Phi2[, i], xout = x)$y/delta)
            if (sum(basis2[, i]) < 0) {
                basis2[, i] <- -basis2[, i]
                coeff2[, i] <- -coeff2[, i]
            }
        }
        colnames(basis2) <- paste("phi", order + (1:m), sep = "")
    }
    return(list(basis = basis, coeff = coeff, varprop = varprop,
        weights = w/mean(w), v = v, basis2 = basis2, coeff2 = coeff2,
        mean.se = axse))
}
