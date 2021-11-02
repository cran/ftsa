# Setting up the functions

# Defining the inner product on L^2([a,b])
inner.product = function(f, g, du) {drop(crossprod(f, g))*du}

# Defining the L2-norm
L2norm = function(f, du) {sqrt(inner.product(f, f, du))}

################
# main function
################

super_fun = function(Y, lag_max, B, alpha, du, p, m, u, select_ncomp, dimension)
{
    n = N = ncol(Y)
    # d0 is the hypothesized dimension
    ## RUN THE BOOTSTRAP IN ORDER TO SELECT d0 PROPERLY
    # d0 = dimension
    
    # Defining Ybar
    Ybar = rowMeans(Y)
    
    # Defining the deviation function Ydev = Y - Ybar, which is used as an input
    # by the function inner.product in constructing the matrix Kstar
    Ydev = Y - Ybar
    
    # 3. Creating the matrix 'Kstar'
    
    # Building the 'core' matrices of Kstar. Below we define the matrices core (n x n),
    # Kstar.core0 [(n-p) x (n-p)] and the array Kstar.core [(n-p) x (n-p) x p].
    # We have that
    #   Kstar = (Kstar.core[,,1] + ... + Kstar.core[,,p])%*%Kstar.core0
    # Where
    #   Kstar.core0[t,s] = <Y[,t],Y[,s]>, t,s in {1,...,n-p}
    #   Kstar.core[t,s,k] = <Y[,t+k],Y[,s+k]>, t,s in {1,...,n-p} and k in {1,...,p}
    # Thus, the matrix core, defined by
    #   core[t,s] = <Y[,t],Y[,s]>, t,s in {1,...,n}
    # contains all the relevant information regarding Kstar.core0 and Kstar.core,
    # and may be used as a building block for the latter matrices.
    core = inner.product(Ydev,Ydev, du)
    Kstar.core0 = core[1:(n-p),1:(n-p)]
    Kstar.core = array(0,c(n-p,n-p,p))
    for (k in 1:p) Kstar.core[,,k] = core[(k+1):(n-(p-k)),(k+1):(n-(p-k))]
    
    # Summing the matrices in 'Kstar.core'
    Kstar.sum = matrix(0,nrow=n-p,ncol=n-p)
    for (k in 1:p) Kstar.sum = Kstar.sum + Kstar.core[,,k]
    
    # Defining Kstar
    Kstar = (n-p)^(-2) * Kstar.sum %*% Kstar.core0
    
    # 4. Eigen-analisys
    
    # Getting the eigenvalues and eigenvectors from 'Kstar'
    
    # Storing the eigenvalues; length(thetahat)=n-p
    thetahat = eigen(Kstar)$values
    
    # Storing the eigenvectors; each column of gammahat corresponds to
    # one eigenvector; dim(gammahat)=[(n-p) x (n-p)]
    gammahat = eigen(Kstar)$vectors
    
    # Defining a tolerance level for 'testing' if the imaginary part of the
    # eigenvalues and eigenvectors is zero
    tol = 10^(-4)
    
    # Checking if there are any complex eigenvalues (among the first eleven)
    for (j in 1:10){
      if (abs(Im(thetahat[j]))>tol) print("Complex eigenvalue found.")
    }
    
    thetahat = Re(thetahat)
    thetahat = sort(thetahat, index.return=TRUE, decreasing=TRUE)
    thetahat.index = thetahat$ix
    thetahat = thetahat$x
    
    # Ordering the eigenvectors accordingly
    gammahat.temp = matrix(0, nrow=nrow(gammahat), ncol=ncol(gammahat))
    for(j in 1:(length(thetahat)))
    {
      gammahat.temp[,j] = gammahat[,thetahat.index[j]]
    }
    gammahat = gammahat.temp
    
    # Storing the original eigenvalues and eigenvectors
    thetahat.old = thetahat
    gammahat.old = gammahat
    
    # Selecting the number of components
    if(select_ncomp == "TRUE")
    {
      bs.pvalues = vector("numeric", lag_max)
      
      # Building the estimated functions Yhat with dimension d0, d0 in {1,2,...,lag_max}.
      # See the main code for an explanation of this section.
      
      for(d0 in 1:lag_max)
      {
        thetahatH0 = Re(thetahat.old[d0+1])
        thetahat = Re(thetahat.old[1:d0])
        gammahat = Re(gammahat.old[,1:d0])
        psihat.root = Ydev[,1:(n-p)] %*% gammahat
        psihat = matrix(0, nrow = m, ncol = d0)
        for(i in 1:d0) psihat[,i] = psihat.root[,i]/L2norm(psihat.root[,i], du = du); rm(i)
        
        etahat = inner.product(psihat, Ydev, du = du)
        Yhat = Ybar + psihat %*% etahat
        
        Yhat.fix = Yhat
        Yhat.fix[Yhat.fix < 0] = 0
        for (t in 1:N) Yhat.fix[,t] = Yhat.fix[,t]/(sum(Yhat.fix[,t])*du); rm(t)
        
        epsilonhat = Y - Yhat.fix
        
        # Let bootstrap begin
        
        sampleCols = function(A)
        {
          idx = sample(1:ncol(A))
          M = matrix(0, nrow = nrow(A), ncol = ncol(A))
          for(i in 1:length(idx))
          {
            M[,i] = A[,idx[i]]
          }
          return(M)
        }
        
        bs.thetahat = vector("numeric", B)
        bs = list()
        for(i in 1:B)
        {
          # This section of the code repeats c-main.r, however inside the list object bs
          # Building the bootstrap 'observed' function Y
          bs$epsilon = sampleCols(epsilonhat)
          bs$Y = Yhat.fix + bs$epsilon
          
          bs$Ybar = rowMeans(bs$Y)
          bs$Ydev = bs$Y - bs$Ybar
          
          bs$core = inner.product(bs$Ydev, bs$Ydev, du = du)
          bs$Kstar.core0 = bs$core[1:(n-p),1:(n-p)]
          bs$Kstar.core = array(0,c(n-p,n-p,p))
          for(k in 1:p) bs$Kstar.core[,,k] = bs$core[(k+1):(n-(p-k)),(k+1):(n-(p-k))]
          bs$Kstar.sum = matrix(0, nrow=n-p, ncol=n-p)
          for(k in 1:p) bs$Kstar.sum = bs$Kstar.sum + bs$Kstar.core[,,k]
          bs$Kstar = (n-p)^(-2) * bs$Kstar.sum %*% bs$Kstar.core0
          bs$thetahat = eigen(bs$Kstar)$values
          bs$thetahat = sort(bs$thetahat, decreasing = TRUE)
          
          # Storing the (d0+1)-th eigenvalue obtained in the i-th bootstrap loop
          bs.thetahat[i] = Re(bs$thetahat[d0+1])
          # print(paste('d0 = ', d0, '; Loop = ', i, sep = ""))
        }
        # Storing the p-values
        bs.pvalues[d0] = sum(bs.thetahat >= thetahatH0)/B
      }        
      d0 = head(which(bs.pvalues < alpha))[1] + 1
      if(is.na(d0)) d0 = 1
    }
    if(select_ncomp == "FALSE")
    {
      d0 = dimension
    }
    # 5. Defining the estimator Yhat
    
    # Storing only the d0 largest eigenvalues and the associated eigenvectors
    thetahat = Re(thetahat.old[1:d0]) # length(thetahat) = d0
    gammahat = Re(gammahat.old[,1:d0]) # dim(gammahat) = [(n-p) x d0]
    
    # Defining the eigenfunctions psihat.root. These are the functions
    # given in equation (2.13)
    #   psihat.root[u,j] = gammahat[1,j]*Ydev[u,1] + ... + gammahat[n-p,j]*Ydev[u,n-p]
    # Note that dim(psihat.root) = [m x d0]. These functions are not
    # necessarily orthonormal.
    psihat.root = Ydev[,1:(n-p)]%*%gammahat
    
    # Normalizing (tested for orthogonality; already orthogonal to each other)
    psihat = matrix(0,nrow=m,ncol=d0)
    #psihat = matrix(0,nrow=length(alpha_grid),ncol=d0)
    for (i in 1:d0) psihat[,i] = psihat.root[,i]/L2norm(psihat.root[,i], du)
    
    # Defining etahat (d0 x n). We have that
    #   etahat[j,t]=<Ydev[,t],psihat[,j]>.
    # This is a d0-vector process. Time varies columnwise.
    etahat = inner.product(psihat,Ydev, du)
    
    # Defining Yhat (m x n). The (u,t)-th element of Yhat is given by
    #   Yhat[u,t] = Ybar[u] + psihat[u,1]*etahat[1,t] + ... + psihat[u,d0]*etahat[d0,t]
    Yhat = Ybar + psihat%*%etahat
    
    # # Note that no restrictions were made upon Yhat, so it may happen that
    # sum(Yhat[,t])*du!=1 and/or Yhat[u,t]<0.
    # We define Yhat.fix to meet these restrictions.
    Yhat.fix = Yhat
    Yhat.fix[Yhat.fix < 0] = 0
    for (t in 1:N) Yhat.fix[,t] = Yhat.fix[,t]/(sum(Yhat.fix[,t])*du); rm(t)
    
    # Defining epsilonhat (m x n)
    epsilonhat = Y - Yhat
    if(select_ncomp == "TRUE")    
    {
      return(list(Y = Y, Ybar = Ybar, thetahat = thetahat, gammahat = gammahat, psihat = psihat,
                  etahat = etahat, Yhat = Yhat, Yhat.fix = Yhat.fix, epsilonhat = epsilonhat, u = u, 
                  d0 = d0, bs.pvalues = bs.pvalues))
    }
    else
    {
      return(list(Y = Y, Ybar = Ybar, thetahat = thetahat, gammahat = gammahat, psihat = psihat,
                  etahat = etahat, Yhat = Yhat, Yhat.fix = Yhat.fix, epsilonhat = epsilonhat, u = u, 
                  d0 = d0))
    }
}
