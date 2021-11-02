ffunopare.knn <-
function(RESPONSES, CURVES, PRED, neighbour,..., kind.of.kernel = "quadratic", 
                          semimetric = "deriv")
{
################################################################
# Performs functional nonparametric prediction (regression) when both
# the response and the predictor are random curves via the functional kernel estimator.
# A global kNN bandwidth is given by the user.
#    "RESPONSES" matrix containing the response curves (row by row)
#    "CURVES" matrix containing the explanatory curves dataset (row by row)
#             used for the estimating stage
#    "PRED" matrix containing new curves stored row by row
#           used for computing predictions
#    "neighbour" number of nearest neighbours used by the estimator
#    "..." arguments needed for the call of the function computing
#          the semi-metric between curves
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
        p1 <- ncol(SEMIMETRIC1)
    n1 <- ncol(SEMIMETRIC1)
        if(neighbour >= n1)
                stop(paste("try a smaller number of neighbour \n than ", neighbour))
        bandwidth.knn1 <- 0
        for(j in 1:p1) {
                Sem <- SEMIMETRIC1[, j]
                knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)]]
                bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
        }
        KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
        KERNEL1[KERNEL1 < 0] <- 0
        KERNEL1[KERNEL1 > 1] <- 0
        diag(KERNEL1) <- 0
        Denom1 <- apply(KERNEL1, 2, sum)
    NUM1 <- t(KERNEL1) %*% RESPONSES
        RESPONSES.estimated <- NUM1/Denom1
    Cv.estimated <- sum((RESPONSES.estimated - RESPONSES)^2)/(n1*p1)
    if(twodatasets) {
        SEMIMETRIC2 <- sm(CURVES, PRED, ...)
        Bandwidth2 <- 0
        p2 <- ncol(SEMIMETRIC2)
                bandwidth.knn2 <- 0
                for(j in 1:p2) {
                        Sem <- SEMIMETRIC2[, j]
                        knn.to.band <- Sem[order(Sem)[neighbour:(neighbour + 1)
                                ]]
                        bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
                }
                KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
                KERNEL2[KERNEL2 < 0] <- 0
                KERNEL2[KERNEL2 > 1] <- 0
                Denom2 <- apply(KERNEL2, 2, sum)
            NUM2 <- t(KERNEL2) %*% RESPONSES
                RESPONSES.predicted <- NUM2/Denom2
        return(list(Estimated.values = RESPONSES.estimated,
            Predicted.values = RESPONSES.predicted, Cv =
            Cv.estimated))
    }else {
        return(list(Estimated.values = RESPONSES.estimated, Cv =
            Cv.estimated))
    }
}
