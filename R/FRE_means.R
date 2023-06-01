FRE_means <-
function(data_pop1,data_pop2,n_year
                    ,n_prefectures,n_age,n_populations=2,row_par){
  Mu_hat<-mu_hat(data_pop1,data_pop2,n_year
                 ,n_prefectures,n_age,n_populations=2,row_par)
  all_data<-array(0,dim=c(n_prefectures,n_populations,n_year,n_age))
  FRE_means<-matrix(0,nrow=n_age,ncol=n_prefectures)
  for (i in 1:n_prefectures) {
    alpha_hat<-rep(0,n_age)
    all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
    all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
    for (j in 1:n_populations) {
      for (k in 1:n_year) {
        alpha_hat=alpha_hat+all_data[i,j,k,]
      }
    }
    FRE_means[,i]<-(1/(n_year*2))*alpha_hat-Mu_hat
  }
  
  return(FRE_means=FRE_means)
  
}
