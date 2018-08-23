kweights <- function (x, kernel = c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", 
    "Quadratic Spectral"), normalize = FALSE) 
{
    kernel <- match.arg(kernel)
    if (normalize) {
        ca <- switch(kernel, Truncated = 2, Bartlett = 2/3, Parzen = 0.539285, 
            `Tukey-Hanning` = 3/4, `Quadratic Spectral` = 1)
    }
    else ca <- 1
    switch(kernel, Truncated = {
        ifelse(ca * x > 1, 0, 1)
    }, Bartlett = {
        ifelse(ca * x > 1, 0, 1 - abs(ca * x))
    }, Parzen = {
        ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5, 1 - 6 * (ca * 
            x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
    }, `Tukey-Hanning` = {
        ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x))/2)
    }, `Quadratic Spectral` = {
        y <- 6 * pi * x/5
        ifelse(x < 1e-04, 1, 3 * (1/y)^2 * (sin(y)/y - cos(y)))
    })
}
