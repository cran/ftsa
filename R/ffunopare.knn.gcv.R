ffunopare.knn.gcv <-
function(RESPONSES, CURVES, PRED, ..., Knearest=NULL, kind.of.kernel = "quadratic", 
                              semimetric = "deriv")
{
    ################################################################
    # Performs functional nonparametric prediction (regression) when both  
    # the response and the predictor are random curves via the functional kernel estimator. 
    # A global bandwidth (i.e. number of neighbours) is selected by a 
    # cross-validation procedure.
    #    "RESPONSES" matrix containing the response curves (row by row) 
    #    "CURVES" matrix containing the explanatory curves dataset (row by row) 
    #             used for the estimating stage
    #    "PRED" matrix containing new curves stored row by row
    #           used for computing predictions
    #    "..." arguments needed for the call of the function computing 
    #          the semi-metric between curves
    #    "Knearest" vector giving the grid of nearest neighbours used in the 
    #               estimator ; if "Knearest"=NULL (default), a grid is 
    #               automatically built
    #    "kind.of.kernel" the kernel function used for computing of 
    #                     the kernel estimator; you can choose 
    #                     "indicator", "triangle" or "quadratic (default)
    #    "semimetric" character string allowing to choose the function 
    #                 computing the semimetric;  you can select 
    #                 "deriv" (default), "fourier", "hshift" and "pca"
    # Returns a list containing:
    #    "Estimated.values" vector containing estimated reponses for 
    #                        each curve of "CURVES"
    #    "Predicted.values" if PRED different from CURVES, this vector 
    #                       contains predicted responses for each 
    #                       curve of PRED
    #    "Bandwidths" vector containing the global data-driven bandwidths  
    #                 for each curve in the matrix "CURVES"
    #    "knearest.opt" optimal number of neighbours
    #    "Cv" cross-validation criterion computed with "CURVES" and "RESPONSES"
    ################################################################
    if(is.vector(RESPONSES)) RESPONSES <- as.matrix(RESPONSES)
    if(is.vector(PRED)) PRED <- as.matrix(t(PRED))
    testfordim <- sum(dim(CURVES)==dim(PRED))==2
    twodatasets <- T
    if(testfordim) twodatasets <- sum(CURVES==PRED)!=prod(dim(CURVES))
    sm <- get(paste("semimetric.", semimetric, sep = ""))
    if(semimetric == "mplsr") stop("semimetric option mlpsr not allowed!")
    SEMIMETRIC1 <- sm(CURVES, CURVES, ...)
    kernel <- get(kind.of.kernel)
    n1 <- ncol(SEMIMETRIC1)
    if(is.null(Knearest))
    {
        step <- ceiling(n1/100)
        if(step == 0) step <- 1
        Knearest <- seq(from = 5, to = n1 %/% 2, by = step)	
        # the vector Knearest contains the sequence of the 
        # k-nearest neighbours used for computing the optimal bandwidth
    }
    kmax <- max(Knearest)	
    p <- ncol(CURVES)
    RESPONSES.estimated <- matrix(0, nrow = n1, ncol = p)
    Bandwidth.opt <- 0
    HAT.RESP <- matrix(0, nrow = n1, ncol = length(Knearest))
    BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
    lKnearest <- length(Knearest)
    HAT.RESP <- array(0,c(n1,lKnearest,p))
    for(i in 1:n1) 
    {
        Norm.diff <- SEMIMETRIC1[, i]	
        # "norm.order" gives the sequence k_1, k_2,... such that
        # dq(X_{k_1},X_i) < dq(X_{k_2},X_i) < ...
        Norm.order <- order(Norm.diff)	
        # "zz" contains dq(X_{k_2},X_i), dq(X_{k_3},X_i),..., 
        # dq(X_{j_{kamx+2}},X_i)
        zz <- sort(Norm.diff)[2:(kmax + 2)]	
        # BANDWIDTH[i, l-1] contains (dq(X_{j_l},X_i) + 
        # dq(X_{j_l},X_i))/2 for l=2,...,kmax+2
        BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
        z <- zz[ - (kmax + 1)]
        ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
        UMAT <- ZMAT/BANDWIDTH[i,  ]
        KNUM <- kernel(UMAT)
        KNUM[col(KNUM) > row(KNUM)] <- 0
        Kdenom <- apply(KNUM[Knearest,  ], 1, sum)
        WEIGHTS <- KNUM[Knearest,  ]/Kdenom
        Ind.curves <- Norm.order[2:(kmax + 1)]
        HAT.RESP[i,,] <- WEIGHTS %*% RESPONSES[Ind.curves,]
    }
    CRITARR <- array(0,c(n1,p,lKnearest))
    for(i in 1:n1)
    {
        CRITARR[i,,] <- (t(HAT.RESP[i,,]) - RESPONSES[i,])^2
    }
    Criterium <- apply(CRITARR, 3, sum)
    index.opt <- order(Criterium)[1]
    RESPONSES.estimated <- HAT.RESP[, index.opt,]
    knearest.opt <- Knearest[index.opt]
    Bandwidth.opt <- BANDWIDTH[, knearest.opt]
    Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p)
    if(twodatasets) 
    {
        SEMIMETRIC2 <- sm(CURVES, PRED, ...)
        Bandwidth2 <- 0
        n2 <- ncol(SEMIMETRIC2)
        for(k in 1:n2) {
          Sm2k <- SEMIMETRIC2[, k]
          Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
        }
        KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
        KERNEL[KERNEL < 0] <- 0
        KERNEL[KERNEL > 1] <- 0
        Denom <- apply(KERNEL, 2, sum)
        NUM <- t(KERNEL) %*% RESPONSES
        RESPONSES.predicted <- NUM/Denom
        return(list(Estimated.values = RESPONSES.estimated, 
                    Predicted.values = RESPONSES.predicted, Bandwidths = 
                      Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                      Cv.estimated))
    }
    else 
    {
        return(list(Estimated.values = RESPONSES.estimated, Bandwidths
                    = Bandwidth.opt, knearest.opt = knearest.opt, Cv = 
                      Cv.estimated))
    }
}
