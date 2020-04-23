farforecast <- function(object, h = 10, var_type = "const", Dmax_value,
                        Pmax_value, level = 80, PI = FALSE)
{
    order_select = method.FPE(object = object, D = Dmax_value, var_type = var_type, Pmax = Pmax_value)

    order = order_select[2]
    ftsm_object = ftsm(y = object, order = order)
    if(requireNamespace("vars", quietly = TRUE))
    {
        var_pred = predict(vars::VAR(ftsm_object$coeff[, 2:(order + 1)], p = order_select[1],
                            type = var_type), n.ahead = h, ci = level/100)
    }
    else
    {
        stop("Please install vars")
    }
    qconf <- qnorm(0.5 + level/200)
    meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = order)
    for(i in 1:order)
    {
        var_fit_pred = var_pred$fcst[[i]]
        meanfcast[, i] = var_fit_pred[, 1]
        varfcast[, i] = ((var_fit_pred[, 3] - var_fit_pred[, 2])/(2 * qconf))^2
    }
    point_fore = ftsm_object$basis[, 2:(order + 1)] %*% t(meanfcast) + ftsm_object$basis[, 1]
    x = as.numeric(rownames(object$y))
    rownames(point_fore) = x
    colnames(point_fore) = 1:h
    point_fore_fts = fts(x, point_fore, yname = "Forecasts", xname = object$xname)
    if (PI == TRUE)
    {
        n.curve = ncol(object$y)
        L = max(round(n.curve/5), order)
        insample_fore = matrix(NA, nrow(object$y), (ncol(object$y) - L))
        for(i in 1:(ncol(object$y) - L))
        {
            dum = ftsm(fts(object$x, object$y[, 1:(L + i - 1)]), order = order)
            dum_coeff = dum$coeff[, 2:(order + 1)]
            var_pred = predict(vars::VAR(dum_coeff, lag.max = nrow(dum_coeff) - 2, type = var_type), n.ahead = 1)
            meanfcast = matrix(NA, nrow = 1, ncol = order)
            for(j in 1:order)
            {
                var_fit_pred = var_pred$fcst[[j]]
                meanfcast[, j] = var_fit_pred[, 1]
            }
            insample_fore[, i] = dum$basis[, 2:(order + 1)] %*% t(meanfcast) + dum$basis[, 1]
        }
        insample_test = object$y[, (L + 1):ncol(object$y)]
        resi = insample_test - insample_fore
        lb_resi = apply(resi, 1, quantile, (100 - level)/200, na.rm = TRUE)
        ub_resi = apply(resi, 1, quantile, (100 + level)/200, na.rm = TRUE)
        lb = point_fore + lb_resi
        ub = point_fore + ub_resi
        colnames(lb) = colnames(ub) = 1:h
        PI_lb = fts(x, lb, yname = "Lower bound", xname = object$xname)
        PI_ub = fts(x, ub, yname = "Upper bound", xname = object$xname)
        return(list(point_fore = point_fore_fts, order_select = order_select,
                    PI_lb = PI_lb, PI_ub = PI_ub))
    }
    else
    {
        return(list(point_fore = point_fore_fts, order_select = order_select))
    }
}
