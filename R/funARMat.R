funARMat <-
function(refinement, eta_star_val)
{
    Mat = funKernel(refinement)
    res <- matrix(nrow = refinement, ncol = ncol(eta_star_val))
    res[,1] = eta_star_val[,1]
    for(ik in 2:(ncol(eta_star_val)))
    {
        res[,ik] = funIntegral(refinement, Mat, eta_star_val[,(ik-1)]) + eta_star_val[,ik]
    }
    return(res)
}
