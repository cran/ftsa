ftsmweightselect <- function(data, ncomp = 6, ntestyear, errorcriterion = c("mae", "mse", "mape"))
{
    errorcriterion = match.arg(errorcriterion)
    time_index = as.numeric(colnames(data$y))
    trainyear  = max(time_index) - ntestyear
    testdat    = extract(data, direction = "time",
                         timeorder = (max(time_index) - ntestyear + 1):max(time_index))$y
    if(errorcriterion == "mae")
    {
        optifunction = function(beta)
        {
            traindata = matrix(NA,length(data$x),ntestyear)
            for(i in 1:ntestyear)
            {
                traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
                                    timeorder = time_index[1]:trainyear), order = ncomp,
                                    weight = TRUE, beta = beta), h = 1)$mean$y
                trainyear = trainyear + 1
            }
            return(mae(traindata, testdat))
        }
        weit = optimize(optifunction,c(0,1))$minimum
    }
    if(errorcriterion == "mse")
    {
        optifunction = function(beta)
        {
            traindata = matrix(NA,length(data$x),ntestyear)
            for(i in 1:ntestyear)
            {
                traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
                                          timeorder = time_index[1]:trainyear), order = ncomp,
                                          weight = TRUE, beta = beta), h = 1)$mean$y
                trainyear = trainyear + 1
            }
            return(mse(traindata, testdat))
        }
        weit = optimize(optifunction,c(0,1))$minimum
    }
    if(errorcriterion == "mape")
    {
        optifunction = function(beta)
        {
              traindata = matrix(NA,length(data$x),ntestyear)
              for(i in 1:ntestyear)
              {
                  traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
                                            timeorder = time_index[1]:trainyear), order = ncomp,
                                            weight = TRUE, beta = beta), h = 1)$mean$y
                  trainyear = trainyear + 1
              }
              return(mape(traindata, testdat))
        }
        weit = optimize(optifunction,c(0,1))$minimum
    }
    return(weit)
}

