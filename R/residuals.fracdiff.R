
# Residuals from arfima() or fracdiff()

residuals.fracdiff <- function(object, ...)
{
    if(!is.null(object$residuals))   # Object produced by arfima()
        return(object$residuals)
    else                             # Object produced by fracdiff()
    {
        if (is.element("x", names(object))) 
            x <- object$x
        else 
            x <- eval.parent(parse(text=as.character(object$call)[2]))
        y <- diffseries(x - mean(x), d=object$d)
        fit <- arima(y, order=c(length(object$ar),0,length(object$ma)), include.mean=FALSE, fixed=c(object$ar,object$ma))
        return(residuals(fit))
    }
}
