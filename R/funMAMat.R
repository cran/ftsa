funMAMat <-
function(refinement, eta_star_val)
{
    Mat = MAphi(refinement, 1.5)
    res <- matrix(nrow = refinement, ncol = ncol(eta_star_val)-1)
    for(ik in 2:(ncol(eta_star_val)))
    {
        res[,ik-1] = funIntegral(refinement, Mat, eta_star_val[,(ik-1)]) + eta_star_val[,ik]
    }
    return(res)
}
