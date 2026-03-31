BrownMat <-
function(N, refinement)
{
    mat <- matrix(nrow = refinement, ncol = N)
    c <- 1
    while(c <= N)
    {
        vec <- BM(N = refinement-1)
        mat[,c] <- vec
        c <- c+1
    }
    return(mat)
}
