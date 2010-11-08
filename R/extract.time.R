extract.time = function (data, timeorder) 
{
    x = data$x
    y = data$y
    index=vector(,length(timeorder))
    if (length(timeorder) == 1){
        index = which(as.numeric(colnames(y)) == timeorder)
    }
    if (length(timeorder) > 1){
        for(i in 1:length(timeorder)){
            index[i] = which(as.numeric(colnames(y)) == timeorder[i])
        }    
    }
    newdata = as.matrix(y[, index])
    ftsdata = fts(x, newdata, xname = data$xname, yname = data$yname)
    colnames(ftsdata$y) = as.character(timeorder)
    return(ftsdata)
}
