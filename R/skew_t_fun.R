skew_t_fun <- function(data, gridpoints, M = 5001)
{
    n = nrow(data)

    # check if input data are densities
    if(all(round(diff(apply(data, 1, sum))) == 0))
    {
        m = gridpoints
    }
    else
    {
        m = seq(min(data), max(data), length.out = M)
    }

    # Parameter estimation

    para_est = sapply(1:n, function(ik) sstdFit(data[ik,])$estimate)
    den_skewed_t = matrix(NA, M, n)
    for(iw in 1:n)
    {
        den_skewed_t[,iw] = dsstd(x = m, mean = para_est[1,iw], sd = para_est[2,iw],
                                  nu = para_est[3,iw], xi = para_est[4,iw])
    }

    # Re-arrange parameters

    para_est_new = matrix(NA, 4, n)

    # log of variances

    para_est_new[2,] = log(para_est[2,]^2)

    # log of transformed df

    x = vector("numeric", length(para_est[3,]))
    for(i in 1:length(para_est[3,]))
    {
        x[i] = pt(q = -2, df = para_est[3,i])
    }
    para_est_new[3,] = log(x)

    # mean and skewness

    para_est_new[1,] = para_est[1,]
    para_est_new[4,] = para_est[4,]

    # differencing log of variance and differencing log of transformed df

    parest_skewt_chen = matrix(nrow = 4, ncol = ncol(para_est) - 1)
    parest_skewt_chen[2,] = diff(para_est_new[2,])
    parest_skewt_chen[3,] = diff(para_est_new[3,])
    parest_skewt_chen[1,] = para_est_new[1,-1]
    parest_skewt_chen[4,] = para_est_new[4,-1]
    rownames(parest_skewt_chen) = c("mean", "diff(log(variance))", "diff(log(transformation))", "skewness")

    # Back-transformation of log of transformed df

    A = exp(para_est_new[3,])
    df_cal = vector()
    for(i in 1:length(A))
    {
        df_cal[i] = uniroot(f = function(df){pt(-2, df) - A[i]},
                            interval = c(min(para_est[3,])-1, max(para_est[3,])+1))$root
    }

    # forecasting four parameters

    VAR_order = VARselect(t(parest_skewt_chen))$selection[1]
    VAR_fore = predict(VAR(t(parest_skewt_chen), p = VAR_order), n.ahead = 1)
    VAR_fore_parest = sapply(1:4, function(ik) VAR_fore$fcst[[ik]][1])

    VAR_fore_parest_new = vector("numeric", 4)
    VAR_fore_parest_new[2] = sqrt(exp(VAR_fore_parest[2] + tail(para_est_new[2,], 1)))

    A_fore = exp(VAR_fore_parest[3] + tail(para_est_new[3,],1))
    dum = try(uniroot(f = function(df){pt(-2, df) - A_fore},
                          interval = c(min(para_est[3,])-1, max(para_est[3,])+1))$root, silent = TRUE)
    if(is(dum, "try-error"))
    {
        print("Unit root problem")
        VAR_fore_parest_new[3] = tail(df_cal,1)
    }
    else
    {
        VAR_fore_parest_new[3] = dum
    }
    VAR_fore_parest_new[1] = VAR_fore_parest[1]
    VAR_fore_parest_new[4] = VAR_fore_parest[4]

    skewed_t_den_fore = dsstd(x = m, mean = VAR_fore_parest_new[1], sd = VAR_fore_parest_new[2],
                              nu = VAR_fore_parest_new[3], xi = VAR_fore_parest_new[4])
    return(list(m = m, skewed_t_den_fore = skewed_t_den_fore))
}
