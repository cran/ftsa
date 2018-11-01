fplsr <- function (data, order = 6, type = c("simpls", "nipals"), unit.weights = TRUE, 
    weight = FALSE, beta = 0.1, interval = FALSE, method = c("delta", 
        "boota"), alpha = 0.05, B = 100, adjust = FALSE, backh = 10) 
{
    type = match.arg(type)
    rawdata = t(data$y)
    n = dim(rawdata)[1]
    Xtrain = rawdata[1:(n - 1), ]
    Ytrain = rawdata[2:n, ]
    Xtest = as.numeric(rawdata[n, ])
    if (interval == FALSE) {
        if (type == "simpls") {
            if (unit.weights == TRUE) 
            {
                output = unitsimpls(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                colnames(fitted) = rownames(Xtrain)
                residuals = t(Ytrain) - fitted
                
                Ypred_mat = as.matrix(output$Ypred)
                colnames(Ypred_mat) = 1:ncol(Ypred_mat)
                
                Xtrain_mat = as.matrix(colMeans(Xtrain))
                colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
                
                Ytrain_mat = as.matrix(colMeans(Ytrain))
                colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                           y1 = as.numeric(colnames(Xtrain)), 
			                     ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			                     y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                           Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
			                     B = output$B, P = output$P, Q = output$Q, T = output$T, R = output$R, 
			                     fitted = fts(1:dim(Xtrain)[2], fitted, xname = data$xname, yname = "Fitted values"), 
                           residuals = fts(1:dim(Xtrain)[2], residuals, xname = data$xname, yname = "Residual"), 
                           meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
                           meanY = fts(1:dim(Ytrain)[2], Ytrain_mat , xname = data$xname, yname = data$yname), 
                           call = match.call())
                return(structure(out, class = "fm"))
            }
            else {
                output = simpls(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
                fitted = t(output$T %*% t(output$Q)) + colMeans(Ytrain)
                colnames(fitted) = rownames(Xtrain)
                residuals = t(Ytrain) - fitted
                
                Ypred_mat = as.matrix(output$Ypred)
                colnames(Ypred_mat) = 1:ncol(Ypred_mat)
                
                Xtrain_mat = as.matrix(colMeans(Xtrain))
                colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
                
                Ytrain_mat = as.matrix(colMeans(Ytrain))
                colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
                out = list(x1 = as.numeric(rownames(Xtrain)), 
                           y1 = as.numeric(colnames(Xtrain)), 
			                     ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
			                     y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
                           Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
			                     B = output$B, P = output$P, Q = output$Q, T = output$T, R = output$R, 
			                     fitted = fts(1:dim(Xtrain)[2], fitted, xname = data$xname, yname = "Fitted values"), 
                           residuals = fts(1:dim(Xtrain)[2], residuals, xname = data$xname, yname = "Residual"), 
                           meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
                           meanY = fts(1:dim(Ytrain)[2], Ytrain_mat, xname = data$xname, yname = data$yname), 
                           call = match.call())
                return(structure(out, class = "fm"))
            }
        }
        else {
            output = nipals(Xtrain, Ytrain, Xtest, order, weight = weight, beta = beta)
            
            Ypred_mat = matrix(output$Ypred, dim(Ytrain)[2], )
            colnames(Ypred_mat) = 1:ncol(Ypred_mat)

            Xtrain_mat = as.matrix(colMeans(Xtrain))            
            colnames(Xtrain_mat) = 1:ncol(Xtrain_mat)
            
            Ytrain_mat = as.matrix(colMeans(Ytrain))
            colnames(Ytrain_mat) = 1:ncol(Ytrain_mat)
            
            fitted_mat_value = t(output$fitted.values[, , order])
            colnames(fitted_mat_value) = rownames(Xtrain)
            
            residual_mat_value = t(output$residuals[, , order])
            colnames(residual_mat_value) = rownames(Xtrain)
            
            out = list(x1 = as.numeric(rownames(Xtrain)), y1 = as.numeric(colnames(Xtrain)), 
      				         ypred = fts(1:dim(Xtrain)[2], t(Xtrain), xname = data$xname, yname = data$yname),
                       y = fts(1:dim(Ytrain)[2], t(Ytrain), xname = data$xname, yname = data$yname), 
      				         Ypred = fts(1:dim(Ytrain)[2], Ypred_mat, xname = data$xname, yname = data$yname), 
      				         P = output$P, Q = output$Q, B = output$B, T = output$T, R = output$R, 
      				         meanX = fts(1:dim(Xtrain)[2], Xtrain_mat, xname = data$xname, yname = data$yname), 
      				         meanY = fts(1:dim(Ytrain)[2], Ytrain_mat, xname = data$xname, yname = data$yname), 
      				         Yscores = output$Yscores, projection = output$projection, 
      				         fitted = fts(1:dim(Xtrain)[2], fitted_mat_value, xname = data$xname, yname = "Fitted values"), 
      				         residuals = fts(1:dim(Xtrain)[2], residual_mat_value, xname = data$xname, yname = "Residual"), 
      				         Xvar = output$Xvar, Xtotvar = output$Xtotvar, 
                       call = match.call())
            return(structure(out, class = "fm"))
        }
    }
    else {
        fplsrPI(t(Xtrain), t(Ytrain), Xtest, order, method = method, 
            alpha = alpha, B = B, weight = weight, beta = beta, 
            adjust = adjust, backh = backh)
    }
}
