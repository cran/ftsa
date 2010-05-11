
# Fitted values from arfima() or fracdiff()

fitted.fracdiff <- function(object, ...)
{
    if(!is.null(object$fitted))      # Object produced by arfima()
        return(object$fitted)
    else
    {
        if (is.element("x", names(object))) 
            x <- object$x
        else 
            x <- eval.parent(parse(text=as.character(object$call)[2]))
        return(x-residuals(object))
    }
}
