extract.x = function (data, xorder) 
{
    x = data$x
    y = data$y
    index=vector(,length(xorder))
    if (length(xorder) == 1){
        index = which(as.numeric(rownames(y)) == xorder)
    }
    if (length(xorder) > 1){
        for(i in 1:length(xorder)){
            index[i] = which(as.numeric(rownames(y)) == xorder[i])
        }    
    }
    newdata = as.matrix(y[index,])
    newx = x[index]
    ftsdata = fts(newx, newdata, xname = data$xname, yname = data$yname)
    rownames(ftsdata$y) = as.character(xorder)
    return(ftsdata)
}
