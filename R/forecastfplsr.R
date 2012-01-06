forecastfplsr <-
function(object, components = 2, h = 20)
{
xname = object$xname
yname = object$yname
foreca = matrix(,length(object$x),h)
for(i in 1:h)
{
fore = fplsr(object, order = components)$Ypred
foreca[,i] = fore$y
newdata = cbind(object$y, fore$y)
colnames(newdata) = as.numeric(colnames(object$y))[1]:(max(as.numeric(colnames(object$y))) + 1)
object = fts(object$x, newdata)
}
forecafts = fts(object$x, foreca, xname = xname, yname = yname)
return(forecafts)
}
