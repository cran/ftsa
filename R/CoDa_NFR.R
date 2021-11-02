CoDa_NFR <-
function(dat, normalize = c("TRUE", "FALSE"), constant)
{
    n_year = nrow(dat)
    n_age = ncol(dat)
    
    # standardize life table death to sum to 1
    
    dat_center = sweep(dat, 1, apply(dat, 1, sum), "/")
    
    # geometric mean
    
    alpha_x = vector("numeric", n_age)
    for(ik in 1:n_age)
    {
        alpha_x[ik] = geometric.mean(dat_center[,ik])
    }
    
    # standardization (closure operation)
    
    f_x_t = matrix(NA, n_year, n_age)
    if(normalize == TRUE)
    {
        for(ik in 1:n_year)
        {
            f_x_t[ik,] = (dat[ik,]/alpha_x)/sum(dat[ik,]/alpha_x)
        }
    }
    if(normalize == FALSE)
    {
        for(ik in 1:n_year)
        {
            f_x_t[ik,] = dat[ik,]/alpha_x
        }
    }
    
    # geometric mean and log-ratio transformation
    
    g_t = vector("numeric", n_year)
    h_x_t = matrix(NA, n_year, n_age)
    if(normalize == TRUE)
    {
        for(ik in 1:n_year)
        {
            g_t[ik] = geometric.mean(f_x_t[ik,])
            h_x_t[ik,] = log(f_x_t[ik,]/g_t[ik])
        }
    }
    if(normalize == FALSE)
    {
        for(ik in 1:n_year)
        {
          h_x_t[ik,] = log(f_x_t[ik,])
        }
    }
  
    # NFR estimation and forecasting
    
    fore_val = as.numeric(ffunopare.knn.gcv(RESPONSES = h_x_t[2:n_year,], CURVES = h_x_t[1:(n_year-1),],
                                 PRED = h_x_t[n_year,], semimetric = "L2")$Predicted.values)
    
    # Inverse clr transformation
    
    f_x_t_star_fore = exp(fore_val)/sum(exp(fore_val))
    d_x_t_star_fore = (f_x_t_star_fore * alpha_x)/sum((f_x_t_star_fore * alpha_x))
    return(d_x_t_star_fore * constant)
}
