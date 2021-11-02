GAEVforecast <-
function(data,q,d.loc.max=10,d.logscale.max=10){

  Tt = dim(data)[1]
  ngrid = dim(data)[2]

  ########################################################################
  # Step 1: select the basis number using leave-one-out cross-validation # 
  ########################################################################
  train = data[1:(Tt-1),]
  test = data[Tt,]
  KLD_test = c()

  for (k1 in 3: d.loc.max){
    for (k2 in 3:d.logscale.max){
      #print(paste('k1=',k1,'k2=',k2))
      para.est = GAEVpara(train,k1,k2)
      para.pred = GAEVfore(para.est)
    
      # the mode of the estimated GAEV model
      gev_mode = function(i){
        mean = para.pred$loc$est[i]
        scale = para.pred$scale$est[i]
        shape = para.pred$shape
        if (shape != 0){
          return (mean+scale*(sign(1+shape)*(abs(1+shape)^(-shape))-1)/shape)
        } else {
          return (mean)
        }
      }
      
      test_mode.pred = apply(as.matrix(1:ngrid),1,gev_mode)
      
      # compare the test data with the MODE of from the estimated GAEV using KLD
      test_KLD = KLD(test_mode.pred,test)$mean.sum.KLD
      KLD_test = rbind(KLD_test,c(k1,k2,test_KLD))

      rm(para.est,para.pred)
    }
  }
  kdf= KLD_test[which.min(KLD_test[,3]),1:2]

  ######################################
  # Step 2: model fitting and forecast #
  ######################################
  para.est = GAEVpara(data,kdf[1],kdf[2])
  para.pred = GAEVfore(para.est)
  
  
  # density forecast
  dens.pred = matrix(0,ncol=ngrid,nrow=length(q))
  if (length(q)==1){
    dens.pred[1,] = apply(as.matrix(1:ngrid),1,
                      function(n){return(dgev(q,loc=para.pred$loc$est[n],scale=para.pred$scale$est[n],shape=para.pred$shape))}) 
  } else{
    for (i in 1:length(q)){
      dens.pred[i,] = apply(as.matrix(1:ngrid),1,
                             function(n){return(dgev(q[i],loc=para.pred$loc$est[n],scale=para.pred$scale$est[n],shape=para.pred$shape))}) 
    }
  }
  rownames(dens.pred) = q
  
  return(list(kdf.location = kdf[1],kdf.logscale = kdf[2],
              basis.location = para.est$loc$basis, basis.logscale = para.est$logscale$basis,
              para.location.pred = para.pred$loc$est, para.scale.pred = para.pred$scale$est, para.shape.pred = para.pred$shape,
              density.pred = dens.pred))
}
