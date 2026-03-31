H_estimate_RS <-
function(fun_dat, band_const = 1, choice_score = c("first", "average"), 
                          est_method = c("classical", "modified"))
{
    est_method = match.arg(est_method)
    T = ncol(fun_dat)
    m = min(T - 1, nrow(fun_dat)) * band_const

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
    
    prop_eigen_value = cumsum(eigen_decomp$values/sum(eigen_decomp$values))
    ncomp = head(which(prop_eigen_value >= 0.95), 1)
   
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
    
    # calculating R_n and S_n values
    
    R_value = var_value = vector("numeric",length(score))
    for(ik in 1:length(score))
    {
        R_value[ik] = sum(score[1:ik] - mean(score))
        var_value[ik] = (score[ik] - mean(score))^2
    }
    R_n = max(R_value) - min(R_value)
    
    S_n_classic = sqrt(sum(var_value)/length(score))
    if(est_method == "classical")
    {
        R_over_S = R_n/S_n_classic
        tilde_H1 = log(R_over_S)/log(length(score))
        d_est = tilde_H1 - 0.5
    }
    if(est_method == "modified")
    {
        rho_hat = acf(score, type = "correlation", plot = FALSE)$acf[2,,1]
        q = floor((length(score) * 3 / 2)^(1/3) * ((2 * rho_hat/(1 - rho_hat^2))^(2/3)))
    
        rho_cov = acf(score, type = "covariance", plot = FALSE)$acf[2:(q+1),,1]
        weight = weight_multi = vector("numeric", q)
        for(ik in 1:q)
        {
            weight[ik] = 1 - ik/(q + 1)
            weight_multi[ik] = weight[ik] * rho_cov[ik] 
        }
        S_n = sqrt(sum(var_value)/length(score) + 2 * sum(weight_multi))
        tilde_H1 = R_n/(S_n * sqrt(length(score))) # for hypothesis test only
        d_est = tilde_H1 - 0.5
    }    
    
    # alpha_estimate
    
    alpha = 1.5 - tilde_H1
    return(list(C = C, score = score, H_est = tilde_H1, d_est = d_est, alpha_est = alpha, 
                ncomp = ncomp, prop_eigen_value = prop_eigen_value, 
                eigen_function = as.numeric(eigen_function)))
}
