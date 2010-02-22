Arima = function (x, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA), xreg = NULL, 
 include.mean = TRUE, include.drift = FALSE, transform.pars = TRUE, fixed = NULL, init = NULL, 
  method = c("CSS-ML", "ML", "CSS"), n.cond, optim.control = list(), kappa = 1e+06, model = NULL)
{
    if (!is.null(model))
        return(arima2(x, model))
    if (include.drift) {
        drift <- 1:length(x)
        xreg <- cbind(xreg, drift)
    }
    if (is.null(xreg))
        tmp <- stats:::arima(x = x, order = order, seasonal = seasonal,
            include.mean = include.mean, transform.pars = transform.pars,
            fixed = fixed, init = init, method = method, n.cond = n.cond,
            optim.control = optim.control, kappa = kappa)
    else tmp <- stats:::arima(x = x, order = order, seasonal = seasonal,
        xreg = xreg, include.mean = include.mean, transform.pars = transform.pars,
        fixed = fixed, init = init, method = method, n.cond = n.cond,
        optim.control = optim.control, kappa = kappa)
    tmp$x <- x
    tmp$series <- deparse(substitute(x))
    if (!is.null(xreg))
        tmp$call$xreg <- xreg
    return(tmp)
}
