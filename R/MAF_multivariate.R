MAF_multivariate <-
function(data, threshold)
{
    n = ncol(data)
    p = nrow(data)
    m = min(n, p)
    Z = scale(t(data), center = TRUE, scale = FALSE)
    S_z = cov(Z)
    svd_1 = eigen(S_z)
    U = svd_1$vectors
    D = svd_1$values
    if(any(D < 0))
    {
        D = replace(D, which(D<0), 0)
    }
    X = Z %*% U %*% diag(D)^(0.5) %*% t(U)
    diff_X = matrix(NA, (n-1), p)
    for(ik in 1:(n-1))
    {
        diff_X[ik,] = X[ik+1,] - X[ik,]
    }
    S_diff = cov(diff_X)
    svd_2 = eigen(S_diff)
    V = svd_2$vectors[,seq(p,1,by=-1)]
    K = rev(svd_2$values)
    W = U %*% diag(D)^(0.5) %*% t(U) %*% V
    Y = Z %*% W

    basis = ginv(Y) %*% Z
    
    basis_norm_const = rowSums(basis)
    basis_norm = apply(basis, 2, "/", basis_norm_const)
    Y_norm = t(apply(Y, 1, "*", basis_norm_const))
    recon_norm = Y_norm %*% basis_norm
    recon_err = Z - recon_norm
    
    # eigenvalues of S^(-1/2)S_{\Delta}S^(-1/2)
    eigen_values = eigen(S_z^(-0.5) %*% S_diff %*% S_z^(-0.5))$values
    V_0 = sum(eigen_values[1:m])
    kmax_star = length(which(eigen_values > V_0/m))
    kmax = min(kmax_star, 0.1 * m)
    ER = GR = vector("numeric", kmax - 1)
    for(k in 1:(kmax - 1))
    {
        ER[k] = eigen_values[k]/eigen_values[k+1]
        V_1 = sum(eigen_values[(k + 1):m])
        u_star_1 = eigen_values[k]/V_1
        V_2 = sum(eigen_values[(k + 2):m])
        u_star_2 = eigen_values[(k + 1)]/V_2
        GR[k] = log(1 + u_star_1)/log(1 + u_star_2)
    }
    k_ER = which.max(ER)
    k_GR = which.max(GR)
    ncomp_eigen_ratio = max(k_ER, k_GR)
    ncomp_threshold = which.min(cumsum(eigen_values) > threshold)
    return(list(MAF = Y_norm, MAF_loading = basis_norm, Z = Z, recon = recon_norm, 
                recon_err = recon_err, ncomp_threshold = ncomp_threshold, ncomp_eigen_ratio = ncomp_eigen_ratio))
}
