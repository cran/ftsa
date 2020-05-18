stop_time_detect <-
function(data, forecasting_method = c("ets", "arima", "rw"))
{
    #################################
    # Forward FTS forecasting errors
    #################################

    p_dim = nrow(data$y)
    n_dim = ncol(data$y)

    ncomp_select_forward = vector("numeric", n_dim - 2)
    ftsm_mat_forward = err_forward = matrix(NA, p_dim, n_dim - 2)
    for(iw in 2:(n_dim - 1))
    {
        train_data = extract(data, timeorder = 1:iw)
        dum = ER_GR(t(train_data$y))
        ncomp_select_forward[iw - 1] = ftsm_ncomp = max(dum$k_ER, dum$k_GR)
        ftsm_mat_forward[,iw - 1] = ftsm_fore = forecast(ftsm(train_data, order = ftsm_ncomp), h = 1,
                                                                       method = forecasting_method)$mean$y
        err_forward[,iw - 1] = data$y[,(iw+1)] - ftsm_fore
        rm(ftsm_fore); rm(train_data); rm(dum)
    }    

    # Begin with 3 curves to estimate

    colnames(err_forward) = 3:n_dim
    ISE_forward_err = ts(colSums(err_forward^2), start = 3, end = n_dim)

    # Breakpoint detection based on ISFE

    ISE_point_forward_strucchange = as.numeric(na.omit(summary(breakpoints(ISE_forward_err ~ 1))$breakdates[1,]))
    ISE_point_forward_ecp = (3:n_dim)[e.divisive(X = matrix(ISE_forward_err, ncol = 1), k = 1)$estimates[2]]
    
    ##################################
    # Backward FTS forecasting errors
    ##################################

    ncomp_select_backward = vector("numeric", n_dim - 2)
    ftsm_mat_backward = err_backward = matrix(NA, p_dim, n_dim - 2)
    for(iw in 2:(n_dim - 1))
    {  
        train_data = extract(data, timeorder = rev(iw:n_dim))
        colnames(train_data$y) = rev(colnames(train_data$y))
        
        dum = ER_GR(t(train_data$y))
        ncomp_select_backward[iw - 1] = ftsm_ncomp = max(dum$k_ER, dum$k_GR)
        ftsm_mat_backward[,iw - 1] = ftsm_fore = forecast(ftsm(train_data, order = ftsm_ncomp), h = 1, 
                                                          method = forecasting_method)$mean$y
        
        err_backward[,iw - 1] = data$y[,iw - 1] - ftsm_fore
        rm(ftsm_fore); rm(train_data); rm(dum)
    }

    # Begin with 3 curves to estimate

    colnames(err_backward) = 1:(n_dim - 2)
    ISE_backward_err = ts(colSums(err_backward^2), start = 1, end = n_dim - 2)

    # Breakpoint detection based on ISFE
    
    ISE_point_backward_strucchange = as.numeric(na.omit(summary(breakpoints(ISE_backward_err ~ 1))$breakdates[1,]))
    ISE_point_backward_ecp = (1:(n_dim - 2))[e.divisive(X = matrix(ISE_backward_err, ncol = 1), k = 1)$estimates[2]]
    
    return(list(break_points_strucchange = c(ISE_point_forward_strucchange, ISE_point_backward_strucchange),
                break_points_ecp = c(ISE_point_forward_ecp, ISE_point_backward_ecp),
                err_forward  = colSums(err_forward^2),
                err_backward = colSums(err_backward^2),
                ncomp_select_forward  = ncomp_select_forward,
                ncomp_select_backward = ncomp_select_backward))
}
