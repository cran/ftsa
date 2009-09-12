`cm.spline` <-
function (x, y = NULL, n = 3 * length(x), method = "fmm", xmin = min(x), 
    xmax = max(x), gulim = 0) 
{
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    nx <- length(x)
    method <- match(method, c("periodic", "natural", "fmm"))
    if (is.na(method)) 
        stop("spline: invalid interpolation method")
    dx <- diff(x)
    if (any(dx < 0)) {
        o <- order(x)
        x <- x[o]
        y <- y[o]
    }
    if (any(diff(y) < 0)) 
        stop("Data are not monotonic")
    if (method == 1 && y[1] != y[nx]) {
        warning("spline: first and last y values differ - using y[1] for both")
        y[nx] <- y[1]
    }
    z <- .C("spline_coef", method = as.integer(method), n = nx, 
        x = x, y = y, b = double(nx), c = double(nx), d = double(nx), 
        e = double(if (method == 1) nx else 0), PACKAGE = "stats")
    z$y <- z$y - z$x * gulim
    z$b <- z$b - gulim
    z <- hyman.filter(z)
    z$y <- z$y + z$x * gulim
    z$b <- z$b + gulim
    z <- spl.coef.conv(z)
    u <- seq(xmin, xmax, length.out = n)
    .C("spline_eval", z$method, nu = length(u), x = u, y = double(n), 
        z$n, z$x, z$y, z$b, z$c, z$d, PACKAGE = "stats")[c("x", 
        "y")]
}

