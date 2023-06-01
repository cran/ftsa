fMBD <-
function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}
