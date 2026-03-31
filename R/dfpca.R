dfpca <- function (x, order, q) 
{
    xs <- scale(x, scale = FALSE)
    n <- nrow(x)
    m <- ncol(x)
    c1 <- array(0, dim = c(n, m, m))
    c2 <- array(0, dim = c(n, m, m))
    for (h in 1:n) {
        for (k in h:n) {
            c1[h, , ] <- c1[h, , ] + xs[k, ] %*% t(xs[k - h + 
                1, ])
            c2[h, , ] <- c2[h, , ] + xs[k - h + 1, ] %*% t(xs[k, 
                ])
        }
        c1[h, , ] <- c1[h, , ]/n
        c2[h, , ] <- c2[h, , ]/n
    }
    f <- array(0, dim = c(m, m))
    for (h in 2:n) {
        f <- f + (wtflat((h/q), 0.5)) * c1[h, , ] + (wtflat((-h/q), 
            0.5)) * c2[h, , ]
    }
    f <- f + c1[1, , ]
    md <- eigen(f, symmetric = TRUE)
    md$values[md$values < 0] <- 0
    score <- xs %*% md$vectors[, 1:order]
    fitted <- t(score %*% t(md$vectors[, 1:order])) + apply(x, 
        2, mean)
    varprop <- cumsum(md$values/sum(md$values))[1:order]
    basis <- cbind(apply(x, 2, mean), md$vectors[, 1:order])
    return(list(coef = score, fitted = fitted, basis = basis, 
        varprop = varprop, order = order))
}