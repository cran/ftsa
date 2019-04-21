FPCA <- function(a)
{
  if (nrow(a)<=1){
    stop("FPCA cannot performed for less than two functions")
  }
  else
  {
    mu<-colMeans(a)
    mean<-matrix(rep(mu,nrow(a)),nrow(a),ncol(a), byrow=TRUE)
    d<-a-mean
    C<- t(d)%*%d/nrow(d)
    wt<-diag(rep(1,ncol(a)))
    eigfn<-sqrt(wt)%*%C%*%sqrt(wt)
    A<-eigen(eigfn)
    eval<-A$values
    evec<-solve(sqrt(wt))%*%A$vectors
    FPCS<-d%*%evec
    eval1<-eval
    eval1[which(A$eval<=0)]=0
    varprop<-rep(0,ncol(a))
    for (np in 1:ncol(a)){
      varprop[np]<-sum(eval1[1:np])/sum(eval1)
    }
    npc.select<-min(which(varprop>=0.9))
    return(list(FPCs=FPCS,eval=eval,evec=evec,mu=mu, mean=mean, npc.select=npc.select))
  }
}
