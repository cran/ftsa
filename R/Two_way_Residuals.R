Two_way_Residuals <-
function(Y,Both,n_prefectures,n_year,n_age){
  Two_FGE=Both$grand_effect
  Two_FRE=Both$row_effect
  Two_FCE=Both$col_effect
  all_male=Y[,1:n_age]
  all_female=Y[,(n_age+1):(2*n_age)]
  
  residuals_b1<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1)) {
    residuals_b1[i,]=all_male[i,]-Two_FGE
    residuals_b2[i,]=all_female[i,]-Two_FGE
  }
  
  #remove the column
  residuals_b1c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1c)) {
    residuals_b1c[i,]=residuals_b1[i,]-Two_FCE[1,]
    residuals_b2c[i,]=residuals_b2[i,]-Two_FCE[2,]
  }
  
  #Remove the row
  residuals_b1r<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2r<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (j in 1:n_prefectures) {
    residuals_b1r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b1c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
    residuals_b2r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b2c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
  }
  
  #Proof Residuals
  data_c1<-matrix(0,nrow=n_prefectures,ncol=n_age)
  data_c2<-matrix(0,nrow=n_prefectures,ncol=n_age)
  for (j in 1:n_prefectures) {
    data_c1[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[1,]
    data_c2[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[2,]
  }
  
  fixed_data_1<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  fixed_data_2<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  for (i in 1:n_prefectures) {
    fixed_data_1[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c1[i,], simplify=FALSE))
    fixed_data_2[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c2[i,], simplify=FALSE))
  }
  FD<-cbind(fixed_data_1,fixed_data_2)
  Recovered_data_1<-fixed_data_1+residuals_b1r
  Recovered_data_2<-fixed_data_2+residuals_b2r
  RD<-cbind(Recovered_data_1,Recovered_data_2)
  
  
  
  #reconstruction proof
  R1=all(round(all_male,4)==round(Recovered_data_1,4))
  R2=all(round(all_female,4)==round(Recovered_data_2,4))
  R<-c(R1,R2)
  
  return(list(residuals1= residuals_b1r,residuals2=residuals_b2r,rd=RD,R=R,Fixed_comp=FD))
  
}
