extract.time = function(data, timeorder)
{
  x = data$x
  y = data$y
  fts(x, y[,timeorder], xname = data$xname, yname = data$yname)
}
