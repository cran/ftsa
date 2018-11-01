centre <- function(x, type)
{
    switch(type, mean = func.mean(t(x)), var = func.var(t(x)),
                 median = depth.FM_fun(t(x))$median,
                 trimmed = depth.FM_fun(t(x))$mtrim)
}
