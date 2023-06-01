wtflat <- function(x, c)
{
    if(-1 < x & x <= -c)
    {
        return(x/(1-c) + 1/(1-c))
    }
    if(-c < x & x < c)
    {
        return(1)
    }
    if(c < x & x < 1)
    {
        return(x/(c - 1) - 1/(c - 1))
    }
    else return(0)
}
