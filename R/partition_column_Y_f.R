partition_column_Y_f <-
function(Y_iter,n_age){
  P_1=Y_iter[,1:n_age]
  P_2=Y_iter[,(n_age+1):(2*n_age)]
  Partition_col=list( P_1,P_2)
  return(Partition_col)
}
