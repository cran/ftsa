H_estimate_LW_NFTS <-
function(fun_dat, band_const, tuning_const, delta_val, m_const,
                       est_method = c("diff_add", "ELW", "ELW_2S"), 
                       choice_score = c("first", "average"), 
                       taper_choice = c("none", "Velasco", "HC"), 
                       ELW_mean_est = c("mean", "init", "weighted", "none"))
{
    T = ncol(fun_dat)
    p_dim = nrow(fun_dat)
    m = min(T - 1, p_dim) * band_const

    X_bar = rowMeans(fun_dat, na.rm = TRUE)
    center_dat = sweep(fun_dat, 1, X_bar)

    # calculating long-run covariance for a given lag

    gamma_l <- function(lag, T)
    {
        gamma_lag_sum = 0
        if(lag >= 0)
        {
            for(ij in 1:(T-lag))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij+lag)])))
            }
        }
        else
        {
            for(ij in 1:(T - abs(lag)))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij + abs(lag)]) %*%
                                t(as.matrix(center_dat[,ij])))
            }
        }
        return(gamma_lag_sum/T)
    }

    # calculating the sum of long-run covariance for all lags (m equals the total number of grid points)

    C = 0
    index = seq(-m, m, by = 1)
    for(k in 1:length(index))
    {
        C = C + (m - abs(index[k])) * gamma_l(lag = index[k], T = T)
    }

    # eigen decomposition

    eigen_decomp = eigen(C)

    # select the number of component based on 95% total amount of variation

    eigen_values = eigen_decomp$values

    m = min(T, p_dim)
    # V_0 = sum(eigen_values[1:m])
    # kmax_star = length(which(eigen_values >= V_0/m))
    # kmax = min(kmax_star, 0.1 * m)

    # ER = vector("numeric", kmax - 1)
    ER = vector("numeric", m)
    for(k in 1:m)
    {
        if(abs(eigen_values[k]/eigen_values[1]) < delta_val)
        {
            ER[k] = 1
        }
        else
        {
            ER[k] = abs(eigen_values[k+1]/eigen_values[1])/abs(eigen_values[k]/eigen_values[1])
        }
    }
    if(all(ER == 1))
    {
        ncomp = 1
    }
    else
    {
        ncomp = which.min(ER) + 1
    }
    
    # first functional principal component

    if(choice_score == "average")
    {
        if(ncomp == 1)
        {
            eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
            score = as.numeric(eigen_function %*% fun_dat)
        }
        if(ncomp > 1)
        {
            eigen_function = eigen_decomp$vectors[,1:ncomp]

            # averaging principal component scores
            score = colMeans(t(eigen_function) %*% fun_dat)
        }
    }
    if(choice_score == "first")
    {
        eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
        score = as.numeric(eigen_function %*% fun_dat)
    }

    ###########################
    # difference + adding back
    # local.W
    ###########################
    
    if(est_method == "diff_add")
    {
        tuning_parameter = floor(1 + length(score)^(m_const))
        d_est = local.W(data = score, m = tuning_parameter, taper = taper_choice)$d
    
        gamma_m = tuning_const/log(tuning_parameter)
        if(d_est < (1 - gamma_m))
        {
            d_est_final = d_est
            ijk = 1
        }
        else
        {
            ijk = 1
            while(d_est >= (1 - gamma_m))
            {
                d_est = local.W(data = diff(score, differences = ijk), m = tuning_parameter, taper = taper_choice)$d
                # print(ijk)
                ijk = ijk + 1
            }
            d_est_final = d_est + (ijk - 1)
        }
    }
    if(est_method == "ELW_2S")
    {
        tuning_parameter = floor(1 + length(score)^0.65)
        d_est_final = ELW2S(data = score, m = tuning_parameter, trend_order = 0, taper = taper_choice)$d
        ijk = 0
    }
    if(est_method == "ELW")
    {
        tuning_parameter = floor(1 + length(score)^0.65)
        d_est_final = ELW(data = score, m = tuning_parameter, mean.est = ELW_mean_est)$d 
        ijk = 0
    }
    tilde_H1 = d_est_final + 0.5
    alpha = 1.5 - tilde_H1
    return(list(C = C, score = score, H_est = tilde_H1, d_est = d_est_final, alpha_est = alpha,
                ncomp = ncomp, ER = ER, diff_order = ijk, eigen_function = as.numeric(eigen_function)))
}
