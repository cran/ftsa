LQDT_FPCA <-
function(data, gridpoints, h_scale = 1, M = 3001,  m = 5001, lag_maximum = 4, no_boot = 1000,
                 alpha_val = 0.1, p = 5, band_choice = c("Silverman", "DPI"),
                 kernel = c("gaussian", "epanechnikov"),
                 forecasting_method = c("uni", "multi"), varprop = 0.85, fmethod, VAR_type)
{
    N = nrow(data)

    # check if input data are densities
    if(all(trunc(diff(apply(data, 1, sum))) == 0))
    {
        Y = t(data)
        u = gridpoints
        du = u[2] - u[1]
    }
    else
    {
        kernel = match.arg(kernel)
        forecasting_method = match.arg(forecasting_method)

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

        # 2. Discretization
        # Evaluation points
        u = seq(from = min(data), to = max(data), length = m)

        # Interval length
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
    }

    # Renormalize Densities to have integral 1
    n = ncol(Y)
    N = length(u)

    any(apply(Y, 2, function(z) any(z < 0))) # make sure no density estimates are negative
    dens = sapply(1:n, function(i) Y[,i]/trapzRcpp(X = u, Y = Y[,i]))

    # Try Forward transformation
    # number of gridpoints for LQD functions - chosen large here so that 0 isn't too close to the boundary of all supports
    lqd = matrix(0, nrow = M, ncol = n)
    c = vector("numeric", n)
    t = seq(0, 1, length.out = M)
    t0 = u[which.min(abs(u))] # closest value to 0
    for(i in 1:n)
    {
        tmp = dens2lqd(dens = dens[,i], dSup = u, lqdSup = t, t0 = t0, verbose = FALSE)
        lqd[,i] = tmp$lqd
        c[i] = tmp$c
    }

    ################################
    # Reconstruction densities
    # Now Try Backward Transform
    # Try cutting off boundary points with large LQD values (effectively setting the density to zero rather than trying to compute it numerically)
    ################################

    cut = res2 = list()
    for(i in 1:n)
    {
        cut[[i]] = c(0, 0)
        cut[[i]][1] = sum(lqd[1:15,i] > 7)
        cut[[i]][2] = sum(lqd[(M-14):M, i] > 7)
        res2[[i]] = lqd2dens(lqd = lqd[,i], lqdSup = t, t0 = t0, c = c[i], cut = cut[[i]], useSplines = FALSE, verbose = FALSE)
    }
    dens2 = sapply(res2, function(r) approx(x = r$dSup, y = r$dens, xout = u, yleft = 0, yright = 0)$y)

    # Assess loss of mass incurred by boundary cutoff
    totalMass = apply(dens2, 2, function(d) trapzRcpp(X = u, Y = d))

    # Numerical comparison of densities

    L2Diff = sapply(1:n, function(i) sqrt(trapzRcpp(X = u, Y = (dens[,i] - dens2[,i])^2))) # L^2 norm difference
    unifDiff = sapply(1:n, function(i){
                                interior = which(dens2[,i] > 0)
                                max(abs(dens[interior,i] - dens2[interior,i]))})
    # Uniform Metric excluding missing boundary values (due to boundary cutoff)

    ######################
    # Forecasting density
    ######################

    c_fore = as.numeric(inv.logit(forecast(auto.arima(logit(c)), h = 1)$mean))

    # cut-off two boundary points

    t_new = t[2:(M-1)]
    lqd_new = lqd[2:(M-1), ]

    if(forecasting_method == "uni")
    {
        ftsm_fitting = ftsm(y = fts(x = t_new, y = lqd_new), order = 10)
        ftsm_ncomp = head(which(cumsum((ftsm_fitting$varprop^2)/sum((ftsm_fitting$varprop^2))) >= varprop),1)
        den_fore = forecast(object = ftsm(y = fts(x = t_new, y = lqd_new), order = ftsm_ncomp), h = 1, method = fmethod)$mean$y
    }
    if(forecasting_method == "multi")
    {
        dt = diff(t)[1]
        foo_out = super_fun(Y = lqd_new, lag_max = lag_maximum, B = no_boot, alpha = alpha_val, du = dt,
                            p = p, m = (M-2), u = t_new, select_ncomp = "TRUE")

        # read outputs

        Ybar_est = foo_out$Ybar
        psihat_est = foo_out$psihat
        etahat_est = matrix(foo_out$etahat, ncol = n)
        selected_d0 = foo_out$d0

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
        den_fore = Ybar_est + psihat_est %*% etahat_pred_val
    }

    # add two useless boundary points back

    den_fore = as.matrix(c(8, den_fore, 8))
    cut = c(sum(den_fore[1:15,1] > 7), sum(den_fore[(M-14):M,1] > 7))

    res_fore = lqd2dens(lqd = den_fore, lqdSup = t, t0 = t0, c = c_fore, cut = cut, useSplines = FALSE, verbose = FALSE)
    dens_fore = approx(x = res_fore$dSup, y = res_fore$dens, xout = u, yleft = 0, yright = 0)$y
    return(list(L2Diff = L2Diff, unifDiff = unifDiff, density_reconstruct = dens2, density_original = dens,
                dens_fore = dens_fore, totalMass = range(totalMass), u = u))
}
