GAEVfore <-
function(para.est){
  nloc = dim(para.est$loc$para)[2]
  nscale = dim(para.est$logscale$para)[2]
  
  loc_pred_para = matrix(0,nrow=nloc)
  logscale_pred_para = matrix(0,nrow=nscale)
  
  allpara = cbind(para.est$loc$para,para.est$logscale$para,para.est$shape)
  colnames(allpara)=1:dim(allpara)[2]
  
  para.forecast = predict(VAR(allpara,ic = "AIC",type='trend'))$fcst
  
  allpara_pred = sapply(para.forecast,'[[',1)
  
  loc_para_pred = allpara_pred[1:nloc]
  logscale_para_pred = allpara_pred[(nloc+1):(nloc+nscale)]
  
  loc_pred = as.vector(para.est$loc$basis%*%loc_para_pred)
  scale_pred = exp(as.vector(para.est$logscale$basis%*%logscale_para_pred))
  shape_pred = allpara_pred[dim(allpara)[2]]
  
  para.pred = list(loc=list(est=loc_pred,para=loc_para_pred),
                   scale=list(est=scale_pred,para=logscale_para_pred),
                   shape=shape_pred)
  
  return(para.pred)
}
