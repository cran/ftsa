LOS <-
function(z, J=NULL, N=NULL, K1=NULL, K2=NULL, phi1=NULL, phi2=NULL, tlength=1){
  resd <- z
  cross.integral <- matrix(0, K1, K2)
  for(i in 1:K1)
    for(j in 1:K2) 
      cross.integral[i,j] <- sum(phi1[,i]* phi2[,j]) *tlength/N
  int1 <- matrix(0, J, K1)
  int2 <- matrix(0, J, K2)
  for(i in 1:J)   {
    for(j in 1:K1)  {
      int1[ i ,j] <- sum( resd[i,] * phi1[,j] ) * tlength /N
    }
    for(j in 1:K2) {
      int2[ i ,j] <- sum( resd[i,] * phi2[,j] ) * tlength /N    
    }
  }
  s1 <- matrix(0, J, K1)
  s2 <- matrix(0, J, K2)
  design.xi <- ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
  resid <- rep(0, K1)
  for(j in 1:J) {
    index <- j
    resid <- resid + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/J
  }
  index.m <- 1:J
  xi.temp <- design.xi %*% resid
  s1[index.m,] <- matrix(rep(xi.temp, each=J), nrow=J)
  for(j in 1:J) {
    index <- j
    s2[index,] <- int2[index,] - drop( t(cross.integral) %*% xi.temp )
  }
  return( list(scores1.out=s1, scores2.out=s2) )
}
