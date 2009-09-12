extract.x = function(data, xorder)
{
  x = data$x
  y = data$y
  newx = x[xorder]
  fts(newx, y[xorder,], xname = data$xname, yname = data$yname)
}

