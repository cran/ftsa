CoDa_FPCA <-
function(data, normalization, h_scale = 1, m = 5001,
                             band_choice = c("Silverman", "DPI"),
                             kernel = c("gaussian", "epanechnikov"),
                             varprop = 0.99, fmethod)
{
    # check if input data are densities

    if(all(trunc(diff(apply(data, 1, sum))) == 0))
    {
        CoDa_mat = t(data)
    }
    else
    {
        band_choice = match.arg(band_choice)
        kernel = match.arg(kernel)
        # Sample size
        N = nrow(data)

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

        # Initialization parameters
        n = N # Number of daily observations

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

        # correcting to ensure integral Y_t du = 1
        for(t in 1:N)
        {
            Y[,t] = Y[,t]/(sum(Y[,t])*du)
        }

        ###########################
        # Dealing with zero values
        ###########################

        return_density_train_trans <- Y
        return_density_train_transformation = return_density_train_trans * (10^6)
        n_1 = ncol(return_density_train_transformation)
        epsilon = sapply(1:n_1, function(X) max(return_density_train_transformation[,X] - round(return_density_train_transformation[,X], 2)))

        CoDa_mat = matrix(NA, m, n_1)
        for(ik in 1:n_1)
        {
            index = which(round(return_density_train_transformation[,ik], 2) == 0)
            CoDa_mat[,ik] = replace(return_density_train_transformation[,ik], index, epsilon[ik])
            CoDa_mat[-index,ik] = return_density_train_transformation[-index,ik] * (1 - (length(index) * epsilon[ik])/(10^6))/(10^6)
        }
    }

    # CoDa

    c = colSums(CoDa_mat)[1]
    dum = CoDa_recon(dat = t(CoDa_mat), normalize = normalization,
                     fore_method = fmethod, fh = 1, varprop = varprop, constant = c)
    return(dum$d_x_t_star_fore)
}
