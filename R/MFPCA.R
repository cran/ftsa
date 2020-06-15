MFPCA <-
function(y,  M=NULL, J=NULL, N=NULL) 
{


  ###     Calculate overall mean function and the deviation 
  ###     from the overall mean to country specific mean functions 
  
  #### Estimate mu ###
  mu<- colMeans(y)
  ### Estimate eta ###
  eta <- matrix(0, M, N)
  for(m in 1:M) {
    eta[ m,  ] <-  colMeans(y[ (1:J)+(m-1)*J,  ]) - mu
  }
  resd <- matrix(0, nrow=M*J, ncol=N) 
  for(m in 1:M) {
    resd[ (1:J)+(m-1)*J,  ] <- t ( t( y[ (1:J)+(m-1)*J,  ] ) - (mu + eta[m,]) ) 
  }
  
  ###    Estimate Rj
  
  Rj<-matrix(0, J, N)
  for (j in 1:J){
    Rj[j,]<-colMeans(resd[(0:(M-1)*J) + j,])
  }
  
  ###    Estimate Uij
  R<-matrix( rep( t( Rj ) , M ) , ncol = ncol(Rj) , byrow = TRUE )
  Uij<-resd - R

  ###     Estimate the covariance functions for time trend and functional pattern 
  G1 <- cov(eta)
  G2 <- long_run_covariance_estimation(t(Rj))
  G3 <- cov(Uij)
  
  
  ###     Estimate eigen values and eigen functions at two levels by calling the 
  ###     eigen function (in R "base" package) on discretized covariance matrices.
  e1 <- eigen(G1)
  e2 <- eigen(G2)
  e3 <- eigen(G3)
  
  ###    transform the eigenvalues of the covariance matrices into eigenvalues
  ###    of original functions
  fpca1.value <- e1$values 
  fpca2.value <- e2$values 
  fpca3.value <- e3$values 
  ###     Keep only non-negative eigenvalues
  fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
  fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
  fpca3.value <- ifelse(fpca3.value>=0, fpca3.value, 0)
 
 
  ###     Calculate the percentage of variance that are explained by the components
  percent1 <- (fpca1.value)/sum(fpca1.value)
  percent2 <- (fpca2.value)/sum(fpca2.value)
  percent3 <- (fpca3.value)/sum(fpca3.value)
  
  ###     Calculate the ratio of first eigenvalue and the rest eigenvalues
  ratio1<- fpca1.value[1]/fpca1.value
  ratio2<- fpca2.value[1]/fpca2.value
  ratio3<- fpca3.value[1]/fpca3.value
  
  ###     Decide the number of components that are kept at level 1 and 2. The general
  ###     rule is to stop at the component where the cumulative percentage of variance 
  ###     explained is greater than 90%.
  
  K1 <-  max(min(which(cumsum(percent1) > 0.9)), min(which(ratio1>sqrt(M)/log10(M)))-1, 2)
  K2 <-  max(min(which(cumsum(percent2) > 0.9)), min(which(ratio2>sqrt(J)/log10(J)))-1, 2)
  K3 <-  max(min(which(cumsum(percent3) > 0.9)), min(which(ratio3>sqrt(M*J)/log10(M*J)))-1, 2)
  
  ###     estimate eigen vectors for discretized covariance matrices and
  ###     transform them into norm one eigenfunctions
  fpca1.vectors <- e1$vectors[, 1:K1]
  fpca2.vectors <- e2$vectors[, 1:K2]
  fpca3.vectors <- e3$vectors[, 1:K3]
  

  
  for(i in 1:K1) {
    v2 <- fpca1.vectors[,i]
    tempsign <- sum(v2)
    fpca1.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
  }
  for(i in 1:K2) {
    v2 <- fpca2.vectors[,i]
    tempsign <- sum(v2)
    fpca2.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
  }
  for(i in 1:K3) {
    v2 <- fpca3.vectors[,i]
    tempsign <- sum(v2)
    fpca3.vectors[,i] <- ifelse(tempsign<0, -1,1) * v2
  }
  
  s1 <- eta%*%fpca1.vectors
  s2 <- Rj%*%fpca2.vectors
  s3 <- Uij%*%fpca3.vectors
  
  

  ###     Return the results from multilevel FPCA as a list
  return( list(K1=K1, K2=K2,K3=K3, lambda1=fpca1.value, lambda2=fpca2.value, lambda3 = fpca3.value, phi1 = fpca1.vectors, phi2=fpca2.vectors, 
               phi3=fpca3.vectors, scores1=s1, scores2=s2, scores3 = s3, mu=mu, eta=t(eta), Rj =Rj, Uij = Uij) )
}
