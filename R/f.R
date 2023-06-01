f <-
function(partition_row_Y,n_year){
  #combine all dataset
  New_data_1<-list()
  New_data_2<-list()
  for (i in 1:length(partition_row_Y)) {
    ND<-partition_row_Y[[i]]
    New_data_1[[i]]=ND[1:n_year,]
    New_data_2[[i]]=ND[(n_year+1):(2*n_year),]
  }
  P_1<-list.rbind(New_data_1)
  P_2<-list.rbind(New_data_2)
  Y_matched<-cbind(P_1,P_2)
  return(Y_matched)
}
