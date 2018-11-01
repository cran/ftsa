dynupdate <- function (data, newdata = NULL, holdoutdata, method = c("ts", 
    "block", "ols", "pls", "ridge"), fmethod = c("arima", "ar", 
    "ets", "ets.na", "rwdrift", "rw"), pcdmethod = c("classical", 
    "M", "rapca"), ngrid = max(1000, ncol(data$y)), order = 6, 
    robust_lambda = 2.33, lambda = 0.01, value = FALSE, interval = FALSE, 
    level = 80, pimethod = c("parametric", "nonparametric"), 
    B = 1000) 
{
    fmethod = match.arg(fmethod)
    pcdmethod = match.arg(pcdmethod)
    pimethod = match.arg(pimethod)
    if (interval == FALSE) {
        ftsmobject = ftsm(data, order = order, ngrid = ngrid, 
            method = pcdmethod, lambda = robust_lambda)
        coef = ftsmobject$coeff
        base = ftsmobject$basis
        p = dim(data$y)[1]
        fore = matrix(NA, (order + 1), 1)
        n2 = length(newdata)
        if (fmethod == "arima") {
            for (i in 1:(order + 1)) {
                fore[i, ] = forecast(auto.arima(coef[, i]), h = 1)$mean
            }
        }
        if (fmethod == "ets") {
            for (i in 1:(order + 1)) {
                fore[i, ] = forecast(ets(coef[, i]), h = 1)$mean
            }
        }
        if (method == "ts") {
            forecasts = (base %*% fore)[(n2 + 1):p, ]
        }
        if (method == "block") {
            updata = data$y[(dim(as.matrix(newdata))[1] + 1):length(as.numeric(data$y))]
            datamatrix = matrix(c(updata, newdata), p, )
            colnames(datamatrix) = colnames(data$y)
            dummy = forecast.ftsm(ftsm(fts(1:p, datamatrix), 
                order = order, ngrid = ngrid, method = pcdmethod), 
                h = 1, method = fmethod, level = level)
            forecasts = dummy$mean$y[1:length(holdoutdata)]
        }
        else {
            base1 = base[1:n2, ]
            base2 = base[(n2 + 1):p, ]
            if (method == "ols") {
                if (length(newdata) == 1) {
                  ols = t(ginv(t(base1) %*% base1) %*% t(base1)) %*% 
                    newdata
                }
                if (length(newdata) > 1) {
                  ols = ginv(t(base1) %*% base1) %*% t(base1) %*% 
                    newdata
                }
                forecasts = base2 %*% ols
            }
            if (method == "pls") {
                I = diag(dim(base)[2])
                if (length(newdata) == 1) {
                  pls = ginv(t(matrix(base1, nrow = 1)) %*% matrix(base1, 
                    nrow = 1) + lambda * I) %*% (t(matrix(base1, 
                    nrow = 1)) %*% newdata + lambda * fore)
                }
                if (length(newdata) > 1) {
                  pls = ginv(t(base1) %*% base1 + lambda * I) %*% 
                    (t(base1) %*% newdata + lambda * fore)
                }
                forecasts = base2 %*% pls
            }
            if (method == "ridge") {
                I = diag(dim(base)[2])
                if (length(newdata) == 1) {
                  ridg = ginv(t(matrix(base1, nrow = 1)) %*% 
                    matrix(base1, nrow = 1) + lambda * I) %*% 
                    (t(matrix(base1, nrow = 1)) %*% newdata)
                }
                if (length(newdata) > 1) {
                  ridg = ginv(t(base1) %*% base1 + lambda * I) %*% 
                    (t(base1) %*% newdata)
                }
                forecasts = base2 %*% ridg
            }
        }
        if (value == TRUE) {
            return(forecasts)
        }
        else {
            errmse = mse(forecasts, holdoutdata)
            errmae = mae(forecasts, holdoutdata)
            errmape = mape(forecasts, holdoutdata)
            return(list(errormse = errmse, errormae = errmae, 
                errormape = errmape))
        }
    }
    else {
        if(method!="pls" & method!="block")
        {
           stop("Interval forecasts can only be updated using the block or pls method.")
        }
        p = nrow(data$y)
        p2 = (length(newdata) + 1):p
        if (method == "pls") {
            if(pimethod != "nonparametric")
            {
                warning("Interval forecast updates from the PLS method can only be obtained via nonparametric bootstrap.")
            }
            output = plsPI(data = data, newdata = newdata, order = order, B = B, 
                           alpha = (100 - level)/100, lambda = lambda)
            plsPI_mat_forecasts = as.matrix(output$forecasts)
            colnames(plsPI_mat_forecasts) = 1
            
            plsPI_mat_bootsamp = as.matrix(output$bootsamp)
            colnames(plsPI_mat_bootsamp) = 1:B
            
            plsPI_mat_low = as.matrix(output$low)
            plsPI_mat_up  = as.matrix(output$up)
            colnames(plsPI_mat_low) = colnames(plsPI_mat_up) = 1
            return(list(forecasts = fts(p2, plsPI_mat_forecasts, xname = data$xname, yname = data$yname), 
                bootsamp = fts(p2, plsPI_mat_bootsamp, xname = data$xname, yname = data$yname), 
                low = fts(p2, plsPI_mat_low, xname = data$xname, yname = data$yname), 
                up = fts(p2, plsPI_mat_up, xname = data$xname, yname = data$yname)))
        }
        if (method == "block") {
            updata = data$y[(dim(as.matrix(newdata))[1] + 1):length(as.numeric(data$y))]
            datamatrix = matrix(c(updata, newdata), p, )
            colnames(datamatrix) = colnames(data$y)
            dummy = forecast.ftsm(ftsm(fts(1:p, datamatrix), 
                order = order, ngrid = ngrid, method = pcdmethod), 
                h = 1, method = fmethod, level = level, pimethod = pimethod, 
                B = B)
            if (pimethod == "parametric") {
                lb_mat = as.matrix(dummy$lower$y[1:length(holdoutdata)])
                ub_mat = as.matrix(dummy$upper$y[1:length(holdoutdata)])
                colnames(lb_mat) = colnames(ub_mat) = 1
                lb = fts(p2, lb_mat, xname = data$xname, yname = data$yname)
                ub = fts(p2, ub_mat, xname = data$xname, yname = data$yname)
                return(list(low = lb, up = ub))
            }
            if (pimethod == "nonparametric") {
                lb = dummy$lower
                ub = dummy$upper
                boot_samp = dummy$bootsamp[1:length(p2), , ]
                return(list(boot_samp = boot_samp, low = lb, up = ub))
            }
        }
    }
}
