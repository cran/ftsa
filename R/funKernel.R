funKernel <-
function(ref)
{
    Mat = matrix(nrow=ref, ncol=ref)
    for(i in 1:ref)
    {
        for(j in 1:ref)
        {
            Mat[i,j] = 0.34 * exp(0.5 * ((i/ref)^2 + (j/ref)^2))
        }
    }
    return(Mat)
}
