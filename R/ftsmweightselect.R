ftsmweightselect <- function(data, ncomp = 6, ntestyear, errorcriterion = c("mae", "mse", "mape"))
{
    errorcriterion = match.arg(errorcriterion)
    trainyear = max(data$time)-ntestyear
    testdat = extract(data, direction = "time", timeorder = (max(data$time)-ntestyear+1):max(data$time))$y
    if(errorcriterion == "mae")
    {
        optifunction = function(beta)
        {
            traindata = matrix(NA,length(data$x),ntestyear)
            for(i in 1:ntestyear)
            {
                traindata[,i] =  forecast(ftsm(extract(data, direction = "time", 
                                    timeorder = data$time[1]:trainyear), order = ncomp,
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
                                          timeorder = data$time[1]:trainyear), order = ncomp,
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
                                            timeorder = data$time[1]:trainyear), order = ncomp,
                                            weight = TRUE, beta = beta), h = 1)$mean$y
                  trainyear = trainyear + 1
              }
              return(mape(traindata, testdat))
        }
        weit = optimize(optifunction,c(0,1))$minimum
    }
    return(weit)
}

