T_stationary <- function(sample, L = 49, J = 500, reps = 1000, seedlen = 50, 
                         quan=c(.90,.95), Ker1 = FALSE, Ker2 = TRUE)
{	
    xrefine = N = ncol(sample)
    refinement = nrow(sample)
    sampling_points = nrow(sample)-1
    
    if(Ker1)
    {
        K=function(x)
        {
            output = min(1, max(1.1-abs(x),0))
            return(output)
        }
    }
    if(Ker2)
    {
        K=function(x)
        {
            output = min(1, max(2-2*abs(x),0))
            return(output)
        }
    }
    
    basis = create.fourier.basis(c(0,1),L)
    time = c(1:(reps/50))*50
    rep = reps
    stataT = c(1:reps)
    ld = length(quan)
    reject95T = reject9T=matrix(0,reps,ld)
    r_95T = r_9T=matrix(0,reps,ld)
    evals = matrix(0,L,reps)
    dres = rp_95 = rp_9 = matrix(0,reps,ld)
    
    Con = .3416
    F = function(x,y)
    {
        return(Con * exp(.5 * (x^2 + y^2)))
    }
    
    Km = matrix(0,refinement,refinement)
    for(j in c(1:refinement))
    {
        for(k in c(1:refinement))
        {
            Km[j,k]=F(j/(refinement),k/(refinement))
        }
    }
    
    trefine=sampling_points+1
    h1=N^.5
    
    for(ll in c(1:reps))
    {
        seed = matrix(0,sampling_points+1,seedlen)
        seederror = matrix(0,sampling_points+1,seedlen)
        for(j in c(1:seedlen))
        {
            seederror[,j] = BBridge(x = 0, y = 0, t0 = 0, T = 1, N = sampling_points)
        }
        seed[,1] = seederror[,1]
        for(j in c(2:seedlen))
        {
            seed[,j]=(Km %*% seed[,j-1]) * (1/(refinement)) + seederror[,j]
        }
        e = matrix(0,(sampling_points+1),N)
        e[,1] = seed[,seedlen]
        
        # sample = matrix(0,sampling_points+1,N)
        # for(ss in c(1:N))
        # {
        #     sample[,ss] = BBridge(x = 0, y = 0, t0 = 0, T = 1, N = sampling_points)
        # }
        
        for(index in 2:N)
        {
            e[,index]=(Km %*% e[,(index-1)])*(1/(refinement))+sample[,index]
        }
        X1_bar = rowMeans(e)
        data1 = e
        mean_subtracted_X1 = data1 - X1_bar
        gamma_hat = list()
        for(i in 0:(N-1))
        {
            temp_matrix = matrix(rep(0, refinement^2), refinement, refinement)
            for(j in (i+1):N)
            {
                temp_matrix = temp_matrix + (mean_subtracted_X1[,j] %*% t(mean_subtracted_X1[,j-i]))
            }
            gamma_hat[i+1]=list(temp_matrix/N)
        }
        cov_sample1 = gamma_hat[[1]]
        for(index in 1:(N-1))
        {
            cov_sample1 = cov_sample1 + K(index/h1) * (gamma_hat[[index+1]] + t(gamma_hat[[index+1]]))
        }
        Z_matrix = cov_sample1
        
        e1 = list()
        for(index in 1:L)
        {
            e1[index] = list(as.matrix(eval.basis(evalarg = (1:refinement)/refinement, basisobj = basis, Lfdobj=0)[,index]))
        }
        eigenvalues1 = (eigen(Z_matrix)$values)/trefine
        D = matrix(0,L,L)
        for(k in 1:L)
        {
            for(ell in 1:L)
            {
                Integrand_matrix = Z_matrix * (e1[[k]] %*% t(e1[[ell]]))
                D[k,ell] = 1/(refinement^2)*sum(Integrand_matrix)
            }
        }
        eigenpairs = eigen(D)
        eigenvectors = eigenpairs$vec
        eigenvalues = eigenpairs$val
        evals[,ll] = eigenvalues
        d = c(1:ld)
        ind = 0
        stoper = 1
        spot = 1
        while(ind==0)
        {
            while((sum(eigenvalues[c(1:spot)])/sum(eigenvalues)) < quan[stoper])
            {
                spot = spot+1
            }
            d[stoper] = spot
            stoper = stoper+1
            if(stoper == (length(d)+1))
            {
                ind = 1
            }
        }
        dres[ll,] = d
        T = c(1:rep)
        lambda = eigenvalues
        p_9T = p_95T = c(1:length(d))
        for(dd in c(1:length(d)))
        {
            for(k in c(1:rep))
            {
                z=rnorm(d[dd]*J)
                tot=0
                for(n in c(1:d[dd]))
                {
                    sum1 = 0
                    sum1 = sum((z[c(((n-1)*d[dd]+1):((n-1)*d[dd]+J))]/(pi*c(1:J)))^2)
                    tot = tot+lambda[n]*sum1
                }
                T[k] = tot
            } 
            fT9 = function(x)
            {
                ecdf(T)(x)-.9
            }
            fT95 = function(x)
            {
                ecdf(T)(x)-.95
            }
            p_9T[dd] = uniroot(fT9, c(mean(T),max(T)),tol=.001,maxiter=1000)$root
            p_95T[dd] = uniroot(fT95, c(mean(T),max(T)),tol=.001,maxiter=1000)$root
            r_9T[ll,dd]=p_9T[dd]
            r_95T[ll,dd]=p_95T[dd]
        }
        Q <- function(x,t)
        {
            data2 = c(1:sampling_points+1)
            data3 = c(1:sampling_points+1)
            for(j in c(1:sampling_points+1))
            {
                data2[j] = sum(data1[j,c(1:floor(N*x))])
                data3[j] = sum(data1[j,c(1:N)])
            }
            return(((1/N)^.5)*(data2[t] - (floor(N*x)/N)*data3[t]))
        }
        int = sum(((1/sqrt(N))*((data1[,c(1:1)])-(1/xrefine)*rowSums(data1)))^2/(xrefine*trefine))
        for(x in c(2:xrefine))
        {
            int = int + sum(((1/sqrt(N))*(rowSums(data1[,c(1:x)])-(x/xrefine)*rowSums(data1)))^2/(xrefine*trefine))
        }
        stataT[ll] = int
        for(dd in c(1:length(d)))
        {
            if(stataT[ll] > p_9T[dd])
            {
                reject9T[ll,dd] = 1
            }
            if(stataT[ll] > p_95T[dd])
            {
                reject95T[ll,dd] = 1
            }
        }
        print(ll)
    }
    return(list(reject9T = reject9T, reject95T = reject95T))
}
