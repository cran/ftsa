mftsc <- function(X, alpha)
{
    fpc.scores <- list()
    npc<-rep(0,length(X))
    for ( j in 1:length(X)){
        fpca.j<-FPCA(X[[j]])
        npc[j]<-fpca.j$npc.select
    }
    npc.select<-max(npc)
    for ( j in 1:length(X)){
        fpca.j<-FPCA(X[[j]])
        fpc.scores[[j]]<- fpca.j$FPCs[,1:npc.select]
    }
    scoredistance<-matrix(rep(0),length(X),length(X))
    for (i in 1:length(X)){
        for (j in 1:length(X))
        scoredistance[i,j]<-tsse(fpc.scores[[i]], fpc.scores[[j]])
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
    cluster.kmeans<-kmeans(distance.score, 2)
    
    member.ks<- cluster.kmeans$cluster
    member.initial<-cbind(c(1:length(X)), member.ks)
    colnames(member.initial)<-c("ID","cluster")
    member_it1<-member.initial
    membership<-list()
    r<-2
    membership[[1]] <- matrix (0,nrow(member.initial),2)
    colnames(membership[[1]])<-c("ID","cluster")
    membership[[2]] <- member_it1
    id.cluster<-sort(unique(member.initial[,2]))
    rand<-rep(0, r)
    while(max(rand)<=alpha){
        member.reclass<-as.vector(rep(0,nrow(member.initial)))
        id.cluster<-sort(unique(member_it1[,2]))
        index<-list()
        for (c in 1:length(id.cluster)){
            index[[id.cluster[c]]]<-as.vector(which(member_it1[,2] == id.cluster[c]))
        }
        for (p in 1:nrow(member.initial)){
            index.out<-index
            index.out[[member_it1[,2][p]]]<-index.out[[member_it1[,2][p]]][-which(index.out[[member_it1[,2][p]]]==p)]
            id.cluster.out<-sort(unique(member_it1[,2][-p]))
            y.pred<-list()
            for (c in 1:length(id.cluster.out)){
                if (length(index.out[[id.cluster.out[c]]])==1){
                    fpca.data<-FPCA(X[[index.out[[id.cluster.out[c]]]]])
                    npc.select<-fpca.data$npc.select
                    mean.data<- matrix(rep(colMeans(X[[p]]), nrow(X[[p]])),nrow(X[[p]]),ncol(X[[p]]), byrow=TRUE)
                    score.out<-(X[[p]]-fpca.data$mean)%*%fpca.data$evec[,1:npc.select]
                    y.pred[[c]] <- score.out%*%t(fpca.data$evec[,1:npc.select]) + mean.data
                }else{
                    data.out<-list()
                    for (d in 1:length(index.out[[id.cluster.out[c]]])){
                        data.out[[d]]<-X[[index.out[[id.cluster.out[c]]][d]]]
                    }
                    M=length(data.out)
                    J=nrow(X[[1]])
                    N=ncol(X[[1]])
                    y<-do.call(rbind, data.out)
                    MFPCA.data<-dmfpca(y, M=M, J=J, N=N, tlength=1)
                    mean.data<- matrix(rep(colMeans (X[[p]]), nrow(X[[p]])),nrow(X[[p]]),ncol(X[[p]]), byrow=TRUE)
                    eta_i<-colMeans(X[[p]])-MFPCA.data$mu
                    res<-X[[p]]-mean.data
                    score.out<-LOS(z=res, J=J, N=N, K1= MFPCA.data$K1, K2=MFPCA.data$K2,
                    phi1=MFPCA.data$phi1, phi2=MFPCA.data$phi2, tlength=1)
                    y.pred[[c]]<- score.out$scores1.out%*%t(MFPCA.data$phi1) +
                    score.out$scores2.out%*%t(MFPCA.data$phi2) + mean.data
                }
            }
            dev<-matrix(0,length(id.cluster.out))
            for(nc in 1:length(id.cluster.out)){
                dev[nc] <- tsse(X[[p]], y.pred[[nc]])
            }
            member.reclass[p] <- id.cluster.out[which(dev == min(dev))]
        }
        member_it1[,2] <- member.reclass
        membership[[r+1]] <- member_it1
        if (length(unique(member.reclass))==1){
            stop("Objects are too similar to cluster")
        }else{
            rand<-rep(NA, r)
            for (t in 1:r){
                rand[t]<-pdfCluster::adj.rand.index(membership[[r+1]][,2], membership[[t]][,2])
            }
            r=r+1
            print(rand)
        }
    }
    return( list(iteration = r-2, membership = membership, member.final = membership[[r]][,2]) )
}
