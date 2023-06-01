get_median <-
function(data){
  if(nrow(data)==2){
    med_curve=colMeans(data)
  }else{
    depth=fMBD(t(data))
    med=which(depth==max(depth))[1]
    med_curve=data[med,]
  }
  return(med_curve)
}
