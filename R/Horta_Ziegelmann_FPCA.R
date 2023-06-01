Horta_Ziegelmann_FPCA <-
function(data, gridpoints, h_scale = 1, p = 5, m = 5001,
                            kernel = c("gaussian", "epanechnikov"),
                            band_choice = c("Silverman", "DPI"),
                            VAR_type = "both", lag_maximum = 6, no_boot = 1000,
                            alpha_val = 0.10, ncomp_select = "TRUE", D_val = 10)
{
    # check if input data are densities

    if(all(trunc(diff(apply(data, 1, sum))) == 0))
    {
        N = nrow(data)
        Y = t(data)
        u = gridpoints
        du = u[2] - u[1]
    }
    else
    {
        # Sample size
        n = N = nrow(data)

        if (!exists('h_scale')) h_scale = 1
        if(band_choice == "Silverman")
        {
            if(kernel == "gaussian")
            {
                h.hat_5m = sapply(1:N, function(t) 1.06*sd(data[t,])*(length(data[t,])^(-(1/5))))
            }
            if(kernel == "epanechnikov")
            {
                h.hat_5m = sapply(1:N, function(t) 2.34*sd(data[t,])*(length(data[t,])^(-(1/5))))
            }
            h.hat_5m = h_scale * h.hat_5m
        }
        if(band_choice == "DPI")
        {
            if(kernel == "gaussian")
            {
                h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "normal"))
            }
            if(kernel == "epanechnikov")
            {
                h.hat_5m = sapply(1:N, function(t) dpik(data[t,], kernel = "epanech"))
            }
            h.hat_5m = h_scale * h.hat_5m
        }

        # defines gridpoints

        u = seq(from = min(data), to = max(data), length = m)
        du = u[2] - u[1]

        # Creating an (m x n) matrix which represents the observed densities. Y[j,t] is the density at date t evaluated at u[j]
        if(kernel == "gaussian")
        {
            Y = sapply(1:N, function(t) density(data[t,], bw = h.hat_5m[t], kernel = 'gaussian', from = min(data), to = max(data), n = m)$y)
        }
        if(kernel == "epanechnikov")
        {
            Y = sapply(1:N, function(t) density(data[t,], bw = h.hat_5m[t], kernel = 'epanechnikov', from = min(data), to = max(data), n = m)$y)
        }

        # correcting to ensure integral Y_t du = 1
        for(t in 1:N)
        {
           Y[,t] = Y[,t]/(sum(Y[,t])*du)
        }
    }

    # main function

    foo_out = super_fun(Y = Y, lag_max = lag_maximum, B = no_boot, alpha = alpha_val, du = du,
                        p = p, m = m, u = u, select_ncomp = ncomp_select, dimension = D_val)

    # read outputs

    Ybar_est = foo_out$Ybar
    psihat_est = foo_out$psihat
    etahat_est = matrix(foo_out$etahat, ncol = N)
    selected_d0 = foo_out$d0
    thetahat_val = foo_out$thetahat
    if(ncomp_select == "TRUE")
    {
        selected_d0_pvalues = foo_out$bs.pvalues
    }
    else
    {
        selected_d0_pvalues = 10^5
    }

    # VAR forecasting

    score_object = t(etahat_est)
    colnames(score_object) = 1:selected_d0

    if(selected_d0 == 1)
    {
        etahat_pred_val = forecast(auto.arima(as.numeric(score_object)), h = 1)$mean
    }
    else
    {
        VAR_mod = VARselect(score_object)
        etahat_pred = predict(VAR(y = score_object, p = min(VAR_mod$selection[3], 3), type = VAR_type), n.ahead = 1)
        etahat_pred_val = as.matrix(sapply(1:selected_d0, function(t) (etahat_pred$fcst[[t]])[1]))
    }
    Yhat.fix_den = den_fore = Ybar_est + psihat_est %*% etahat_pred_val

    # adjustment

    Yhat.fix_den[Yhat.fix_den < 0] = 0
    Yhat.fix_den = Yhat.fix_den/(sum(Yhat.fix_den)*du)
    return(list(Yhat.fix_den = Yhat.fix_den, u = u, du = du, Ybar_est = Ybar_est,
                psihat_est = psihat_est, etahat_est = etahat_est, etahat_pred_val = etahat_pred_val,
                selected_d0 = selected_d0, selected_d0_pvalues = selected_d0_pvalues, thetahat_val = thetahat_val))
}
