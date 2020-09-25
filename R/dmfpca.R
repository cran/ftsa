dmfpca <- function(y,  M=NULL, J=NULL, N=NULL, tstart=0, tlength=1)
{
    ###     Generate the set of grid points that are equally spaced.
    t <- seq(tstart, tstart + tlength, length=N)
    
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
    
    
    ###     Estimate the three covariance functions: overall covariance G, 
    ###     between covariance Gb and within covariance Gw
    # Gb <- long_run_covariance_estimation(t(Rj))$BT_FT_fix_C0
    Gb <- long_run_covariance_estimation(t(Rj))
    Gw <- cov(Uij)
    
    
    ###     Estimate eigen values and eigen functions at two levels by calling the 
    ###     eigen function (in R "base" package) on discretized covariance matrices.
    e1 <- eigen(Gb)
    e2 <- eigen(Gw)
    
    ###    transform the eigenvalues of the covariance matrices into eigenvalues
    ###    of original functions
    fpca1.value <- e1$values * tlength / N
    fpca2.value <- e2$values * tlength / N 
    ###     Keep only non-negative eigenvalues
    fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
    fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
    ###     Calculate the percentage of variance that are explained by the components
    percent1 <- (fpca1.value)/sum(fpca1.value)
    percent2 <- (fpca2.value)/sum(fpca2.value)
    ###     Calculate the ratio of first eigenvalue and the rest eigenvalues
    ratio1<- fpca1.value[1]/fpca1.value
    ratio2<- fpca2.value[1]/fpca2.value
    ###     Decide the number of components that are kept at level 1 and 2. The general
    ###     rule is to stop at the component where the cumulative percentage of variance 
    ###     explained is greater than 90%.
    
    K1 <-  max(min(which(cumsum(percent1) > 0.9)), min(which(ratio1>sqrt(J)/log10(J)))-1, 2)
    K2 <-  max(min(which(cumsum(percent2) > 0.9)), min(which(ratio2>sqrt(M*J)/log10(M*J)))-1, 2)
    
    ###     estimate eigen vectors for discretized covariance matrices and
    ###     transform them into norm one eigenfunctions
    fpca1.vectors <- e1$vectors[, 1:K1]*sqrt(N/tlength)
    fpca2.vectors <- e2$vectors[, 1:K2]*sqrt(N/tlength)
    
    ###     The eigenfunctions are unique only up to a change of signs.
    ###     Select the signs of eigenfunctions so that the integration over the domain 
    ###     is non-negative
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
    
    #############################################################################
    ###     The following part will estimate the principal component scores
    #############################################################################
    
    ###     First, calculate the inner product (the cosine of the angles) between 
    ###     level 1 eigenfunctions and level 2 eigenfunctions
    cross.integral <- matrix(0, K1, K2)
    for(i in 1:K1)
        for(j in 1:K2) 
            cross.integral[i,j] <- sum(fpca1.vectors[,i]* fpca2.vectors[,j]) *tlength/N
    
    ### Create the coefficient matrix ###
    I.matrix<-diag(K1)
    Z<-cbind(I.matrix, cross.integral)
    
    ###     Next, calculate the inner product of each centered function with the 
    ###     level 1 or level 2 eigenfunctions
    int1 <- matrix(0, M*J, K1)
    int2 <- matrix(0, M*J, K2)
    for(i in 1:(M*J))   {
        for(j in 1:K1)  {
            int1[ i ,j] <- sum( resd[i,] * fpca1.vectors[,j] ) * tlength /N
        }
        for(j in 1:K2) {
            int2[ i ,j] <- sum( resd[i,] * fpca2.vectors[,j] ) * tlength /N    
        }
    }
    
    
    ###     Finally, calculate the principal component scores based on the formulas
    ###     given in the paper.
    
    beta.hat <- MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% t(int1)
    s1.raw<-t(beta.hat[1:K1,])
    s1<- matrix(0, J, K1)
    for (j in 1:J){
        s1[j,]<-colMeans(s1.raw[(0:(M-1)*J) + j,])
    }
    s2<-t(beta.hat[-(1:K1),])
    
    
    
    ###     Return the results from multilevel FPCA as a list
    return( list(K1=K1, K2=K2, lambda1=fpca1.value, lambda2=fpca2.value, phi1=fpca1.vectors, 
                 phi2=fpca2.vectors, scores1=s1, scores2=s2, mu=mu, eta=t(eta)) )
}
