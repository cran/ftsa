ftsmiterativeforecasts <-
function(object, components, iteration = 20)
{
xname = object$xname
yname = object$yname
foreca = matrix(, length(object$x), iteration)
for(i in 1:iteration)
{
fore = forecast(ftsm(object, order = components), h = 1)$mean$y
foreca[,i] = fore
newdata = cbind(object$y, fore)
colnames(newdata) = as.numeric(colnames(object$y))[1]:(max(as.numeric(colnames(object$y))) + 1)
object = fts(object$x, newdata)
}
forecafts = fts(object$x, foreca, xname = xname, yname = yname)
return(forecafts)
}
