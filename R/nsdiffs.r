nsdiffs = function (x, m = frequency(x))
{
    if (m <= 1)
        stop("Non seasonal data")
    chstat <- SD.test(x, m)
    crit.values <- c(0.4617146, 0.7479655, 1.0007818, 1.237535,
        1.462524, 1.69202, 1.9043096, 2.1169602, 2.3268562, 2.5406922,
        2.7391007)
    if (m <= 12)
        D <- as.numeric(chstat > crit.values[m - 1])
    else if (m == 24)
        D <- as.numeric(chstat > 5.098624)
    else if (m == 52)
        D <- as.numeric(chstat > 10.341416)
    else if (m == 365)
        D <- as.numeric(chstat > 65.44445)
    else D <- as.numeric(chstat > 0.269 * m^(0.928))
    return(D = D)
}
