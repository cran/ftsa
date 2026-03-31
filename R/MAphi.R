MAphi <-
function(ref, k)
{
    Mat = matrix(nrow = ref, ncol = ref)
    for(i in 1:ref)
    {
        for(j in 1:ref)
        {
            Mat[i,j] = min(i,j)/ref
        }
    }
    return(k * Mat)
}
