as_flag <- function(x, name)
{
    if (is.character(x))
        x = match.arg(x, c("TRUE", "FALSE"))
    out = as.logical(x)[1L]
    if (is.na(out))
        stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
    out
}