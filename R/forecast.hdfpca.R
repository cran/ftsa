forecast.hdfpca <- function(object, h = 3, level = 80, B = 50, ...)
{
    order = object$order
    r = object$r
    m = object$m
    p = object$p
    n = dim(object$y[[1]])[2]
    resid<- list()
    for(im in 1: m)
    {
        resid[[im]] <- object$y[[im]]-object$fitted[[im]]
    }

    mod.fore = load.fore <- list()
    for(io in 1:order)
    {
        load.fore[[io]] <- array(NA, dim = c(h, r))
        mod.fore[[io]] <- list()
        for(ir in 1:r)
        {
            mod <- auto.arima(object$model$model2[[io]]$coef[, ir]) # forecast the factor loadings
            mod.fore[[io]][[ir]] <- forecast(mod, h)
            load.fore[[io]][,ir] <- mod.fore[[io]][[ir]]$mean
        }
    }
    score.fore <- list()
    for(io in 1:order) # the forecast scores
    {
        score.fore[[io]] <- object$model$model2[[io]]$basis[, 1] + object$model$model2[[io]]$basis[, 2:(1+r)] %*% t(load.fore[[io]])
    }
    score <- list()
    for(im in 1:m) # combine the forecast socres for each population
    {
        score[[im]] <- sapply(score.fore, '[', ((im-1)*h+1):(im*h))
    }
    # reconstruct forecast functions
    fun.fore <- list()
    if(h == 1)
    {
        for(im in 1: m)
        {
            fun.fore[[im]] <- object$model$model1[[im]]$basis[, 1] + object$model$model1[[im]]$basis[, 2:(1+order)] %*% score[[im]]
        }
    }
    else if (h >= (n/2))
    {
        stop('forecast horizon is too big considering the sample size')
    }
    else
    {
        for(im in 1: m)
        {
            fun.fore[[im]] <- object$model$model1[[im]]$basis[, 1] + object$model$model1[[im]]$basis[, 2:(1+order)] %*% t(score[[im]])
        }
    }
    boot <- list()
    for(io in 1:order)
    {
        boot[[io]] <- array(NA, dim = c(r, h, B))
        for(ir in 1: r)
        {
            boot[[io]][ir,,] <- t(sco.resamp(mod.fore[[io]][[ir]], h = h, B = B))
        }
    }
    fsid = boot.beta = boot.curve <- list()
    for(im in 1:m)
    {
        boot.beta[[im]] <- array(NA, dim = c(order, h, B))
        boot.curve[[im]] <- array(NA, dim = c(p, h, B))
    }
    modhd <- object$model$model2
    for(ib in 1: B)
    {
        fsid <- list()
        for(im in 1: m)
        {
            fsid[[im]] <- array(NA, dim = c(p, h, B))
            abc <- list() # bootstrapped score forecast
            for(io in 1:order)
            {
                abc[[io]] <- (modhd[[io]]$basis[, 1] + modhd[[io]]$basis[,2:(r+1)]%*%boot[[io]][,,ib])[im,]
            }
            boot.beta[[im]][,,ib] <- t(sapply(abc, '[', 1:h)) # combine bootstrapped score forecast
            for(ih in 1:h)
            {
                for(ip in 1: p) # resample function residual
                {
                    fsid[[im]][ip, ih, ] <- sample(resid[[im]][ip,], B, replace = T)
                }
            }
            boot.curve[[im]][,,ib] <- object$model$model1[[im]]$basis[,1] + object$model$model1[[im]]$basis[,2:(order+1)]%*%boot.beta[[im]][,,ib] + fsid[[im]][,,ib]
        }
    }
    upper = lower = list()
    for(im in 1:m)
    {
        upper[[im]] <- apply(boot.curve[[im]], c(1, 2), quantile, c((1-level/100)/2, 1-(1-level/100)/2))[1,,]
        lower[[im]] <- apply(boot.curve[[im]], c(1, 2), quantile, c((1-level/100)/2, 1-(1-level/100)/2))[2,,]
    }
    return(structure(list(forecast = fun.fore, upper = upper, lower = lower), class = "forecast.hdfpca"))
}
