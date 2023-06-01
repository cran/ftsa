DFPCA <- function(y)
{
  mu <- colMeans(y)
  sub_mean <- matrix(rep(mu,nrow(y)),nrow(y),ncol(y), byrow=TRUE)
  resd <- y-sub_mean
  G <- long_run_covariance_estimation(t(resd))
  e1 <- eigen(G)
  fpca.value <- e1$values
  fpca.value <- ifelse(fpca.value>=0, fpca.value, 0)
  percent <- (fpca.value)/sum(fpca.value)
  ratio <- fpca.value[1]/fpca.value
  K <-  max(min(which(cumsum(percent) > 0.9)), min(which(ratio>sqrt(nrow(y))/log10(nrow(y))))-1, 2)
  fpca.vectors <- e1$vectors
  FPCS <- resd %*% fpca.vectors
  return(list(score = FPCS, lambda = fpca.value, phi = fpca.vectors,
              mu = mu, npc.select = K, mean = sub_mean))
}
