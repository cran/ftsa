partition_row_Y_f <-
function(Y_iter,column_partition_index,row_partition_index){
  partition_Y_1 <- lapply(1:length(column_partition_index), function(k) {
    Y_iter[,column_partition_index[[k]] ]})
  partition_by_row=list()
  for (h in 1:length(partition_Y_1)) {
    YY=partition_Y_1[[h]]
    partition_by_row[[h]] <- lapply(1:length(row_partition_index), function(k) {
      YY[row_partition_index[[k]], ]})
  }
  Final_partition_by_row<-list()
  for (k in 1:length(partition_by_row[[1]])) {
    PBR_1=partition_by_row[[1]]
    PBR_2=partition_by_row[[2]]
    Final_partition_by_row[[k]]=rbind(PBR_1[[k]],PBR_2[[k]])
  }
  return(Final_partition_by_row)
}
