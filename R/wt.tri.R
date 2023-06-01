wt.tri <- function(x)
{
    if(abs(x)>1)
    {
        return(0)
    }
    else
    {
        return(1 - abs(x))
    }
}
