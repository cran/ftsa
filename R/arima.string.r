arima.string = function (object)
{
    order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
    result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3],
        ")", sep = "")
    if (order[7] > 1 & sum(order[4:6]) > 0)
        result <- paste(result, "(", order[4], ",", order[5],
            ",", order[6], ")[", order[7], "]", sep = "")
    if (is.element("constant", names(object$coef)) | is.element("intercept",
        names(object$coef)))
        result <- paste(result, "with non-zero mean")
    else if (is.element("drift", names(object$coef)))
        result <- paste(result, "with drift")
    else if (order[2] == 0 & order[5] == 0)
        result <- paste(result, "with zero mean")
    else result <- paste(result, "")
    return(result)
}
