hdfpca <- function(y, order, q = sqrt(dim(y[[1]])[2]), r)
{
    mod = score = basis = varprop <- list() # all fpca models
    m <- length(y) # number of populations
    n <- dim(y[[1]])[2] # sample size
    p <- dim(y[[1]])[1] # number of discrete points
    
    for(im in 1:m) # for each population, do dynamic fpca
    {
        mod[[im]] <- dfpca(t(y[[im]]), order, q)
        score[[im]] <- mod[[im]]$coef
        basis[[im]] <- mod[[im]]$basis
        varprop[[im]] <- mod[[im]]$varprop
    }
    
    score.cb <- list()
    for(io in 1:order) # extract the first few pca scores from all populations
    {
        score.cb[[io]] <- sapply(score, '[', ((io-1)*n+1):(io*n))
    }
    
    ## High-dim pca
    modhd = load = factor = fitted = varprophd <- list()
    for(io in 1:order) # further do pca on the combined scores
    {
        modhd[[io]] <- hdpca(score.cb[[io]], r)
        load[[io]] <- modhd[[io]]$coef
        factor[[io]]<- modhd[[io]]$basis
        fitted[[io]] <- t(modhd[[io]]$fitted)
        varprophd[[io]] <- modhd[[io]]$varprop
    }
    score.fit <- list() # extract and combine the fitted scores
    for(im in 1:m)
    {
        score.fit[[im]] <- sapply(fitted, '[', ((im-1)*n+1):(im*n))
    }
    # function reconstruction
    fun.fit <- list()
    for(im in 1: m)
    {
        fun.fit[[im]] <- basis[[im]][, 1] + basis[[im]][, 2:(1+order)] %*% t(score.fit[[im]])
    }
    out = list(y = y, p = p, fitted = fun.fit, m = m, model = list(model1 = mod, model2 = modhd),
               order = order, r = r)
    return(structure(out, class = "hdfpca"))
}
