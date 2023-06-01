mu_hat <- function(data_pop1, data_pop2, n_year, n_prefectures, n_age, n_populations = 2, row_par)
{
  all_data <- array(0, dim = c(n_prefectures, n_populations, n_year, n_age))
  Mu_hat <- rep(0,n_age)
  for (i in 1:n_prefectures)
  {
      all_data[i,1,,]=t(data_pop1[,(n_year*i-(n_year-1)):(n_year*i)])
      all_data[i,2,,]=t(data_pop2[,(n_year*i-(n_year-1)):(n_year*i)])
      for (j in 1:n_populations)
      {
          for (k in 1:n_year)
          {
              Mu_hat=Mu_hat+all_data[i,j,k,]
          }
      }
  }
  final_GE<-1/(n_year*n_populations*n_prefectures)*Mu_hat
  return(final_GE)
}
