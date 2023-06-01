FCE_means <-
function(data_pop1,data_pop2,n_year
                    ,n_prefectures,n_age,n_populations=2,row_par){
  Mu_hat<-mu_hat(data_pop1,data_pop2,n_year
                 ,n_prefectures,n_age,n_populations=2,row_par)
  all_data<-array(0,dim=c(n_prefectures,n_populations,n_year,n_age))
  FCE_means<-matrix(0,nrow=n_age,ncol=n_populations)
  for (j in 1:n_populations) {
    beta_hat<-rep(0,n_age)
    for (i in 1:n_prefectures) {
      all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
      all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
      for (k in 1:n_year) {
        beta_hat=beta_hat+all_data[i,j,k,]
      }
    }
    FCE_means[,j]<-(1/(n_year*n_prefectures))*beta_hat-Mu_hat
  }
  return(FCE_means=FCE_means)
}
