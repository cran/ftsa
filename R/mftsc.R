mftsc<-function(X, alpha){
  dfpc.scores <- list()
  npc<-rep(0,length(X))
  for ( j in 1:length(X)){
    dfpca.dat<-DFPCA(X[[j]])
    npc[j]<-dfpca.dat$npc.select
    dfpc.scores[[j]]<-dfpca.dat$score
  }
  npc.select<-max(npc)
  score.select<-list()
  for ( j in 1:length(X)){
    score.select[[j]]<- dfpc.scores[[j]][,1:npc.select]
  }
  scoredistance<-matrix(rep(0),length(X),length(X))
  for (i in 1:length(X)){
    for (j in 1:length(X))
      scoredistance[i,j]<-tsse(score.select[[i]], score.select[[j]])
  }
  
  distance.score<-as.dist(scoredistance)
  distortion<-rep(0,7)
  for (c in 1:7){
    distortion[[c]]<-kmeans(distance.score, c)$tot.withinss/c
  }
  trans_distort<-rep(0,6)
  for (c in 2:7){
    trans_distort[[c-1]]<-distortion[[c]]^(-1/2)-distortion[[c-1]]^(-1/2)
  }
  trans_distort
  ncluster = which(trans_distort == max(trans_distort))+1
  cluster.kmeans<-kmeans(distance.score, ncluster)
  member.ks<- cluster.kmeans$cluster
  #member.initial<-data.frame(V1=c(1:20),V2=c(rep(1,10),rep(2,10)),V3=as.vector(member.ks))
  member.initial <- data.frame(V1 = c(1:length(X)),
                                 V2 = rep(NA, length(X)),
                                 V3 = as.vector(member.ks))
  colnames(member.initial)<-c("ID","ture membership","cluster")
  
  member_it1<-member.initial
  membership<-list()
  r<-2
  membership[[1]] <- matrix(0.1, length(X), 3)
  colnames(membership[[1]])<-c("ID","ture membership","cluster")
  membership[[2]] <- member_it1
  id.cluster<-sort(unique(member.initial[,2]))
  rand<-rep(0, r)
  while(max(rand)<=alpha){
    member.reclass<-as.vector(rep(0,nrow(member.initial)))
    for (p in 1:nrow(member.initial)){
      Ypred<-list()
      for (k in 1:ncluster){
        #### if k is the same cluster as pth object
        if (member_it1$cluster[p]==k){
          if (length(which(member_it1$cluster==k))==1){
            dfpca.pred<-DFPCA(X[[p]])
            mean.pred<- matrix(dfpca.pred$mu, nrow = nrow(X[[1]]) ,ncol = length(dfpca.pred$mu), byrow=TRUE)
            score.pred<-(X[[p]]-dfpca.pred$mean)%*%dfpca.pred$phi[,1:dfpca.pred$npc.select]
            Ypred[[k]] <- score.pred%*%t(dfpca.pred$phi[,1:dfpca.pred$npc.select]) + mean.pred
          }else{
            ###########Leave one#########
            id.out<-which(member_it1$cluster==k)[-which(which(member_it1$cluster==k)==p)]
            if (length(id.out)==1){
              dfpca.pred<-DFPCA(X[[id.out]])
              mean.pred<- matrix(dfpca.pred$mu, nrow = nrow(X[[1]]) ,ncol = length(dfpca.pred$mu), byrow=TRUE)
              score.pred<-(X[[p]]-dfpca.pred$mean)%*%dfpca.pred$phi[,1:dfpca.pred$npc.select]
              Ypred[[k]] <- score.pred%*%t(dfpca.pred$phi[,1:dfpca.pred$npc.select]) + mean.pred
            }else{
              dat.out<-list()
              for (l in 1:length(id.out)){
                dat.out[[l]]<-X[[id.out[l]]]
              }
              M=length(dat.out)
              J=nrow(dat.out[[1]])
              N=ncol(dat.out[[1]])
              y<-do.call(rbind, dat.out)
              MFPCA.pred<-MFPCA(y, M=M, J=J, N=N)
              grandmean<-matrix(rep(MFPCA.pred$mu, nrow(X[[1]])),
                                ncol = ncol(X[[1]]), byrow = TRUE)
              mean.p<-matrix(rep( colMeans(X[[p]]), nrow(X[[p]])),
                             ncol = ncol(X[[p]]), byrow = TRUE)
              eta.p <- colMeans(X[[p]])-MFPCA.pred$mu
              score1.out <- eta.p%*%MFPCA.pred$phi1
              res<- X[[p]]-matrix(rep(colMeans(X[[p]]), nrow(X[[p]])),nrow(X[[p]]),ncol(X[[p]]), byrow=TRUE)
              score.out<-LOS(z=res, J=J, N=N, K1= MFPCA.pred$K2, K2=MFPCA.pred$K3,
                             phi1=MFPCA.pred$phi2, phi2=MFPCA.pred$phi3)
              score2.out<-score.out$scores1.out
              score3.out<-score.out$scores2.out
              Ypred[[k]] <- grandmean + matrix(rep( score1.out%*%t(MFPCA.pred$phi1), nrow(X[[p]])),
                                               ncol = ncol(X[[p]]), byrow = TRUE) + score2.out%*%t(MFPCA.pred$phi2) + score3.out%*%t(MFPCA.pred$phi3)
            }
          }
          ##### if k is not the same as pth object
        }else{
          if (length(which(member_it1$cluster==k))<=1){
            if (length(which(member_it1$cluster==k))==0){
              Ypred[[k]]<-matrix(0, nrow = nrow(X[[p]]), ncol = ncol(X[[p]]))
            }else{
              dfpca.pred<-DFPCA(X[[which(member_it1$cluster==k)]])
              mean.pred<- matrix(dfpca.pred$mu, nrow = nrow(X[[1]]) ,ncol = length(dfpca.pred$mu), byrow=TRUE)
              score.pred<-(X[[p]]-dfpca.pred$mean)%*%dfpca.pred$phi[,1:dfpca.pred$npc.select]
              Ypred[[k]] <- score.pred%*%t(dfpca.pred$phi[,1:dfpca.pred$npc.select]) + mean.pred
            }
          }else{
            ####MFPCA on this cluster
            id.out<-which(member_it1$cluster==k)
            dat.out<-list()
            for (l in 1:length(id.out)){
              dat.out[[l]]<-X[[id.out[l]]]
            }
            M=length(dat.out)
            J=nrow(dat.out[[1]])
            N=ncol(dat.out[[1]])
            y<-do.call(rbind, dat.out)
            MFPCA.pred<-MFPCA(y, M=M, J=J, N=N)
            grandmean<-matrix(rep(MFPCA.pred$mu, nrow(X[[1]])),
                              ncol = ncol(X[[1]]), byrow = TRUE)
            mean.p<-matrix(rep( colMeans(X[[p]]), nrow(X[[p]])),
                           ncol = ncol(X[[p]]), byrow = TRUE)
            eta.p <- colMeans(X[[p]])-MFPCA.pred$mu
            score1.out <- eta.p%*%MFPCA.pred$phi1
            res<- X[[p]]-matrix(rep(colMeans(X[[p]]), nrow(X[[p]])),nrow(X[[p]]),ncol(X[[p]]), byrow=TRUE)
            score.out<-LOS(z=res, J=J, N=N, K1= MFPCA.pred$K2, K2=MFPCA.pred$K3,
                           phi1=MFPCA.pred$phi2, phi2=MFPCA.pred$phi3)
            score2.out<-score.out$scores1.out
            score3.out<-score.out$scores2.out
            Ypred[[k]] <- grandmean + matrix(rep( score1.out%*%t(MFPCA.pred$phi1), nrow(X[[p]])),
                                             ncol = ncol(X[[p]]), byrow = TRUE) + score2.out%*%t(MFPCA.pred$phi2) + score3.out%*%t(MFPCA.pred$phi3)
          }
        }
      }
      diff<-rep(0, ncluster)
      for (i in 1:ncluster){
        diff[i]<-tsse(as.matrix(Ypred[[i]]), as.matrix(X[[p]]))
      }
      member.reclass[p] <- which.min(diff)
    }
    member_it1[,3] <- member.reclass
    membership[[r+1]] <- member_it1
    if(length(unique(member.reclass))==1){
      rand<-c(rand,1)
    }else{
      rand<-rep(NA, r)
      for (t in 1:r){
        rand[t]<-adj.rand.index(membership[[r+1]][,3], membership[[t]][,3])
      }
    }
    r <- r+1
    print(rand)
  }
  return( list(iteration = r-2, membership = membership, member.final = membership[[r]][,2]) )
}


