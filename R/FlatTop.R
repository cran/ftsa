FlatTop <- function(x)
{
    out <- 0
    if(x < 0.5)
    {
        out <- 1
    }
    else if(x >= 0.5 && x <= 1)
    {
        out <- 2 - 2*x
    }
    else 
    {
        out <- 0
    }
    return(out)
}
