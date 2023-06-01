list.do <-
function (.data, fun, ...) 
{
  do.call(what = fun, args = as.list(.data), ...)
}
