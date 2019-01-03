facf <-
function(fun_data, lag_value_range = seq(0, 20, by = 1))
{
    center_dat = scale(fun_data, center = TRUE, scale = FALSE)
    T = nrow(center_dat)
    gamma_l <- function(lag, T)
    {
        gamma_lag_sum = 0
        if(lag >= 0)
        {
            for(ij in 1:(T-lag))
            {
                gamma_lag_sum = gamma_lag_sum + as.matrix(center_dat[ij,]) %*% t(as.matrix(center_dat[(ij+lag),]))
            }
        }
        else
        {
            for(ij in 1:(T+lag))
            {
                gamma_lag_sum = gamma_lag_sum + as.matrix(center_dat[ij-lag,]) %*% t(as.matrix(center_dat[ij,]))
            }
        }
        return(gamma_lag_sum/T)
    }
    
    rho_val = vector("numeric", length(lag_value_range))
    for(ik in 1:length(lag_value_range))
    {
        lag_value = lag_value_range[ik] 
        if(lag_value == 0)
        {
            rho_val[ik] = NA
        }
        else
        {
            rho_val[ik] = sqrt(sum((gamma_l(lag = lag_value, T = T))^2))/sum(diag(gamma_l(lag = 0, T = T)))
        }
    }
    return(rho_val)
}
