dynupdateweightselect <-
function(data, method  = c("pls", "ridge"), interval = c(0, 10^4), p, backh = 1, errortype = c("mse", "mae", "mape"))
{
  if(missing(method)){
     method = "pls"
  }
  if(missing(errortype)){
     errortype = "mse"
  }
  han = function(lambda){
        n = dim(data$y)[2]
        dimen = dim(data$y)[1]
        errors = vector(, backh)
        if (method == "pls"){
            if(errortype == "mse"){            
              for(j in 1:backh){
                  errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                    holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "pls", fmethod = "ets", error = "mse", lambda = lambda)$error
              }
            }
            if(errortype == "mae"){
               for(j in 1:backh){
                   errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                    holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "pls", fmethod = "ets", error = "mae", lambda = lambda)$error
               }
            }
            if(errortype == "mape"){
               for(j in 1:backh){
                   errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                    holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "pls", fmethod = "ets", error = "mape", lambda = lambda)$error
               }
            }      
        }
        else{
            if(errortype == "mse"){
               for(j in 1:backh){
                   errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                     holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "ridge", fmethod = "ets", error = "mse", lambda = lambda)$error
               }
            }
            if(errortype == "mae"){
               for(j in 1:backh){
                   errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                     holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "ridge", fmethod = "ets", error = "mae", lambda = lambda)$error
               }
            }
            if(errortype == "mape"){
               for(j in 1:backh){
                   errors[j] = dynupdate(data = fts(1:dimen, data$y[, 1:(n-j)]), newdata = data$y[1:p, (n-j+1)],
                                     holdoutdata = data$y[(p+1):dimen, (n-j+1)], method = "ridge", fmethod = "ets", error = "mape", lambda = lambda)$error
               }
            }             
        }
        return(mean(errors))
  }
  optimize(han, interval)
}

