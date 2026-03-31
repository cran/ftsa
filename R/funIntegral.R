funIntegral <-
function(ref, Mat, X)
{	
    Mat = Mat %*% X
    return(Mat/ref)
}
