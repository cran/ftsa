newmodel = function (p, d, q, P, D, Q, constant, results)
{
    n <- nrow(results)
    for (i in 1:n) {
        if (identical(c(p, d, q, P, D, Q, constant), results[i,
            1:7]))
            return(FALSE)
    }
    return(TRUE)
}
