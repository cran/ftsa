Substraction <-
function(data){
  med<-get_median(data)
  new_data=matrix(0,nrow = dim(data)[1],ncol = dim(data)[2])
  for (i in 1:dim(data)[1]) {
    new_data[i,]=data[i,]-med
  }
  return(new_data)
}
