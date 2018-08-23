plot.ftsm <- function (x, components, components.start = 0, xlab1 = x$y$xname, ylab1 = "Basis function",
xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean",
level.lab = "Level", main.title = "Main effects", interaction.title = "Interaction",
basiscol = 1, coeffcol = 1, outlier.col = 2, outlier.pch = 19,
outlier.cex = 0.5, ...)
{
    oldpar <- par(no.readonly = TRUE)
    mean <- is.element("mean", colnames(x$basis))
    level <- is.element("level", colnames(x$basis))
    if(!(mean | level))
    components.start = max(1, components.start)
    m <- mean + level
    order <- ncol(x$basis) - m
    if (missing(components))
    components <- order
    n <- components -  components.start + 1
    par(mfcol = c(2, n))
    if(components.start == 0){
        if (mean)
        plot(x$y$x, x$basis[, "mean"], type = "l", lty = 1, xlab = xlab1,
        ylab = mean.lab, main = main.title, col = basiscol,
        ...)
        if (level)
        plot.ts(x$coeff[, "level"], xlab = xlab2, ylab = level.lab,
        main = ifelse(mean, "", main.title), col = coeffcol,
        ...)
        if (m == 1)
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n",
        xlab = "", ylab = "")
    }
    for (i in max(1, components.start):components) {
        yl1 <- ifelse(n > 1, paste(ylab1, i), ylab1)
        yl2 <- ifelse(n > 1, paste(ylab2, i), ylab2)
        plot(x$y$x, x$basis[, m + i], type = "l", lty = 1,
        xlab = xlab1, ylab = yl1, col = basiscol, ...)
        if (i == 1)
        title(interaction.title)
        plot.ts(x$coeff[, m + i], xlab = xlab2, ylab = yl2,
        col = coeffcol, ...)
        if (!is.null(x$wt)) {
            if (sum(x$wt < 0.1)/length(x$wt) < 0.2) {
                points(time(x$coeff)[x$wt < 0.1], x$coeff[x$wt <
                0.1, i + m], pch = outlier.pch, col = outlier.col,
                cex = outlier.cex)
            }
        }
    }
    par(oldpar)
}

