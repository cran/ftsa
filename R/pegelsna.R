pegelsna <- function (x, upper = c(alpha = 0.3, beta = 0.2, phi = 0.99),
                         lower = c(alpha = 0.01, beta = 0.01, phi = 0.9), model = "AZN")
{
    MSE1 <- function(alpha, x) {
        if (alpha > upper[1])
        return(1e+09)
        if (alpha < 0.01)
        return(1e+09)
        fit <- Arima(x, order = c(0, 1, 1), fixed = alpha - 1)
        return(fit$sigma2)
    }
    MSE2 <- function(alpha, x) {
        if (any(alpha > upper[1:2]))
        return(1e+09)
        if (any(alpha < lower[1:2]))
        return(1e+09)
        theta1 <- alpha[1] + alpha[1] * alpha[2] - 2
        theta2 <- 1 - alpha[1]
        fit <- Arima(x, order = c(0, 2, 2), fixed = c(theta1,
        theta2))
        return(fit$sigma2)
    }
    MSE3 <- function(alpha, x) {
        if (any(alpha > upper))
        return(1e+09)
        if (any(alpha < lower))
        return(1e+09)
        if (alpha[3] < alpha[2])
        return(1e+09)
        theta1 <- alpha[1] + alpha[1] * alpha[2] - 1 - alpha[3]
        theta2 <- (1 - alpha[1]) * alpha[3]
        phi1 <- alpha[3]
        fit <- Arima(x, order = c(1, 1, 2), fixed = c(phi1, theta1,
        theta2))
        return(fit$sigma2)
    }
    
    n <- length(x)
    
    if (model == "AZN") {
        fit1 <- nlm(MSE1, (upper[1] + lower[1]) / 2, x = x)
        fit2 <- nlm(MSE2, (upper[1:2] + lower[1:2]) / 2, x = x)
        fit3 <- nlm(MSE3, (upper + lower) / 2, x = x)
        n <- length(x)
        aic <- c(n * log(fit1$minimum) + 2, n * log(fit2$minimum) +
        4, n * log(fit3$minimum) + 6)
        best <- c("ANN", "AAN", "ADN")[aic == max(aic)]
    }
    else
    best <- model
    
    if (best == "ANN") {
        fit1 <- nlm(MSE1, (upper[1] + lower[1]) / 2, x = x, stepmax = 10)
        fitarima <- Arima(x, order = c(0, 1, 1), fixed = fit1$estimate -
        1)
        fitpar <- c(alpha = fit1$estimate,
        beta = 0,
        gamma = 0,
        phi = 1)
        method <- "Robust SES"
    }
    else if (best == "AAN") {
        fit2 <- nlm(MSE2, (upper[1:2] + lower[1:2]) / 2, x = x)
        theta1 <- fit2$estimate[1] + fit2$estimate[1] * fit2$estimate[2] -
        2
        theta2 <- 1 - fit2$estimate[1]
        fitarima <- Arima(x, order = c(0, 2, 2), fixed = c(theta1,
        theta2))
        fitpar <- c(alpha = fit2$estimate[1],
        beta = fit2$estimate[2],
        gamma = 0,
        phi = 1)
        method <- "Robust Holt's"
    }
    else if (best == "ADN") {
        fit3 <- nlm(MSE3, (upper + lower) / 2, x = x)
        theta1 <- fit3$estimate[1] + fit3$estimate[1] * fit3$estimate[2] -
        1 - fit3$estimate[3]
        theta2 <- (1 - fit3$estimate[1]) * fit3$estimate[3]
        phi1 <- fit3$estimate[3]
        fitarima <- Arima(x, order = c(1, 1, 2), fixed = c(phi1,
        theta1, theta2))
        fitpar <- c(alpha = fit3$estimate[1],
        beta = fit3$estimate[2],
        gamma = 0,
        phi = fit3$estimate[3])
        method <- "Robust Damped Holt's"
    }
    else stop("Unknown model")
    
    fitarima$par <- fitpar
    fitarima$method <- method
    return(fitarima)
}

