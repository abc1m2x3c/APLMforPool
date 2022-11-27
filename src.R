#this is the source function where estimator has a sample weight original.weight
library(maxLik)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("PLM.cpp")

fit.aplm<-function(Z,U,X,original.Weight,PoolID,U.ps,X.pred,U.pred,kernel_type,max.iter=10,plot=FALSE,plot.pred=FALSE){
  Z=as.matrix(Z)
  U=as.matrix(U)
  X=as.matrix(X)
  
  p=ncol(X)
  q=ncol(U)
  N=nrow(Z)
  J=length(unique(PoolID))
  
  #sort all data by PoolID
  index=order(PoolID)
  Z=Z[index]
  U=U[index,]
  X=X[index,]
  original.Weight=original.Weight[index]
  PoolID=PoolID[index]
  
  #identify the start and end location of each pool
  unique_PoolID=unique(PoolID)
  start_ind=rep(0,length(unique_PoolID))
  end_ind=rep(0,length(unique_PoolID))
  count=0
  for (id in unique_PoolID){
    count=count+1
    index=which(PoolID==id)
    start_ind[count]=index[1]
    end_ind[count]=index[length(index)]
  }
  
  Weight=rep(NA,N)
  for (i in 1:length(start_ind)){
    index=start_ind[i]:end_ind[i]
    Weight[index]=original.Weight[index]/sum(original.Weight[index])
  }
  
  #pre calculate some variables to speed up
  #pre calculation for linear
  # H=solve(t(X)%*%X)%*%t(X)
  X.pool=matrix(NA,nrow=J,ncol=p)
  Z.pool=rep(NA,J)
  original.Weight.pool=rep(NA,J)
  for (j in 1:J){
    index=start_ind[j]:end_ind[j]
    X.pool[j,]=colSums(matrix(X[index,]*Weight[index],nrow=length(index),ncol=p))
    Z.pool[j]=Z[index][1]
    original.Weight.pool[j]=sum(original.Weight[index])
  }
  H.pool=solve(t(X.pool*original.Weight.pool)%*%X.pool)%*%t(X.pool*original.Weight.pool)
  
  #pre calculation for nonlinear
  U.pool=matrix(NA,ncol=q,nrow=J)
  for (j in 1:J){
    index=start_ind[j]:end_ind[j]
    for (i in 1:q){
      U.pool[j,i]=sum(U[index,i]*Weight[index])
    }
  }
  
  # #pre define vectors to record results
  record.beta=matrix(NA,nrow=p,ncol=max.iter)
  record.f=array(0,c(N,q,max.iter))
  record.h=matrix(NA,nrow=q,ncol=max.iter)
  
  #initialize
  mu.hat=mean(Z)
  f.hat=matrix(0,nrow=N,ncol=q)
  Y.hat.pre=rep(mu.hat,N)
  beta.hat.pre=rep(0,p)
  cond=TRUE
  iter=0
  
  #algorithm starts
  while (iter==0 | cond==TRUE){
    iter=iter+1
    Y.hat=matrix(rep(NA,N),nrow=N,ncol=1)
    for (j in 1:J){
      index=start_ind[j]:end_ind[j]
      for (i in 1:length(index)){
        Y.hat[index[i]]=(Z[index[i]]-sum(Weight[index[-i]]*Y.hat.pre[index[-i]]))/Weight[index[i]]
      }
    }
    
    #update linear part
    f.hat.pool=matrix(NA,nrow=J,ncol=q)
    for (j in 1:J){
      index=start_ind[j]:end_ind[j]
      f.hat.pool[j,]=colSums(matrix(f.hat[index,]*Weight[index],nrow=length(index),ncol=q))
    } 
    beta.hat=H.pool%*%(Z.pool-rowSums(f.hat.pool))
    X.beta.hat=X%*%beta.hat

    
    #calculate X.beta.hat.pool in preparation for homogeneous pooling in nonlinear part
    X.beta.hat.pool=rep(NA,J)
    for (j in 1:J){
      index=start_ind[j]:end_ind[j]
      X.beta.hat.pool[j]=sum(X.beta.hat[index]*Weight[index])  
    }
    
    #update nonlinear part
    for (i in 1:q){
      if (U.ps[i]=='nonhomo'){ #non homo pooling Ver.
        cv_h<-function(h){
          Y.temp=Y.hat-X.beta.hat-rowSums(as.matrix(f.hat[,-i]))
          res=rep(NA,N)
          for (j in 1:J){
            index=start_ind[j]:end_ind[j]
            Sc.temp=LLS_W_cpp(U[-index,i],h,U[index,i],original.Weight[-index],kernel_type) #S1
            res[index]=(Y.temp[index]-Sc.temp%*%Y.temp[-index])^2
          }
          return (mean(res))
        }
        tt=proc.time()
        print(paste("The ",iter,'th iteration, optimizing h',i,'.',sep=''))
        h=optimize(cv_h,lower=diff(range(U[,i]))/5,upper=diff(range(U[,i]))/2)$minimum
        print(paste('Done. h',i,'=',h,' Time elapse is ',round((proc.time()-tt)[3],2),sep=''))
        S1=LLS_W_cpp(U[,i],h,U[,i],original.Weight,kernel_type)
        S1c=S1-matrix(colMeans(S1),byrow=TRUE,nrow=N,ncol=N)
        f.hat[,i]=S1c%*%(Y.hat-X.beta.hat-rowSums(as.matrix(f.hat[,-i])))
      }else{   #homo pooling Ver.
        f.hat.pool=matrix(NA,nrow=J,ncol=q)
        for (j in 1:J){
          index=start_ind[j]:end_ind[j]
          f.hat.pool[j,-i]=colSums(matrix(f.hat[index,-i]*Weight[index],nrow=length(index),ncol=q-1))
        }
        cv_h<-function(h){
          Y.temp=Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool
          res=rep(NA,J)
          for (j in 1:J){
            Sc.temp=LLS_W_cpp(U.pool[-j,i],h,U.pool[j,i],original.Weight.pool[-j],kernel_type) 
            res[j]=(Y.temp[j]-Sc.temp%*%Y.temp[-j])^2
          }
          return (mean(res))
        }
        tt=proc.time()
        print(paste("The ",iter,'th iteration, optimizing h',i,'.',sep=''))
        h=optimize(cv_h,lower=diff(range(U[,i]))/5,upper=diff(range(U[,i]))/2)$minimum
        print(paste('Done. h',i,'=',h,' Time elapse is ',round((proc.time()-tt)[3],2),sep=''))
        S1=LLS_W_cpp(U.pool[,i],h,U.pool[,i],original.Weight.pool,kernel_type)
        S1c=S1-matrix(colMeans(S1),byrow=TRUE,nrow=J,ncol=J)
        f.hat.pool[,i]=S1c%*%(Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool)  
        for (j in 1:J){
          index=start_ind[j]:end_ind[j]
          f.hat[index,i]=f.hat.pool[j,i]
        }
      }
      record.h[i,iter]=h
    }
    
    Y.hat.pre=X.beta.hat+rowSums(f.hat)
    cond= (iter<max.iter)
    beta.hat.pre=beta.hat
    
    
    #record the result of each iteration
    record.beta[,iter]=beta.hat
    record.f[,,iter]=f.hat
    
    #plot part
    if (plot==TRUE){
      par(mfrow=c(1,q))
      sqerror=0
      for (i in 1:q){
        xlab.name=paste('z',i,sep='')
        ylab.name=paste('f',i,'(z',i,')',sep='')
        if (U.ps[i]=='homo'){
          Y.plot.temp=rep(NA,J)
          U.plot.temp=rep(NA,J)
          for (j in 1:J){
            index=start_ind[j]:end_ind[j]
            Y.plot.temp[j]=sum(Weight[index]*f.hat[index,i])
            U.plot.temp[j]=sum(Weight[index]*U[index,i])
          }
          order.temp=order(U.plot.temp)
          plot(U.plot.temp[order.temp],Y.plot.temp[order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
        }else{
          order.temp=order(U[,i])
          plot(U[,i][order.temp],f.hat[,i][order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
        }
      }
    }
    Sys.sleep(1)
  }
  
  record.beta=as.matrix(record.beta[,1:iter])
  record.f=as.matrix(record.f[,,1:iter])
  record.h=as.matrix(record.h[,1:iter])
  
  #prediction part
  f.pred=matrix(NA,nrow=nrow(as.matrix(U.pred)),ncol=q)
  for (i in 1:q){
    if (U.ps[i]=='nonhomo'){
      S1.pred=LLS_W_cpp(U[,i],record.h[i,iter],U.pred[,i],original.Weight,kernel_type)
      S1.original=LLS_W_cpp(U[,i],record.h[i,iter],U[,i],original.Weight,kernel_type)
      S1c=S1.pred-matrix(colMeans(S1.original),byrow=TRUE,nrow=nrow(f.pred),ncol=N)
      f.pred[,i]=S1c%*%(Y.hat.pre-X.beta.hat-rowSums(as.matrix(f.hat[,-i]))) #use Y.hat.pre
    }
    if (U.ps[i]=='homo'){
      #predict based on homo pooled results
      S1.pred=LLS_W_cpp(U.pool[,i],record.h[i,iter],U.pred[,i],original.Weight.pool,kernel_type)
      S1.original=LLS_W_cpp(U.pool[,i],record.h[i,iter],U.pool[,i],original.Weight.pool,kernel_type)
      S1c=S1.pred-matrix(colMeans(S1.original),byrow=TRUE,nrow=nrow(f.pred),ncol=J)
      f.pred[,i]=S1c%*%(Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool)
    }
  }
  #Y.pred
  Y.pred=X.pred%*%beta.hat+rowSums(as.matrix(f.pred))
  if (plot.pred==TRUE){
    par(mfrow=c(1,q))
    sqerror=0
    for (i in 1:q){
      xlab.name=paste('z',i,sep='')
      ylab.name=paste('f',i,'(z',i,')',sep='')
      if (U.ps[i]=='homo'){
        Y.plot.temp=(f.pred[,i])
        U.plot.temp=(U.pred[,i])
        
        order.temp=order(U.plot.temp)
        plot(U.plot.temp[order.temp],Y.plot.temp[order.temp],type='l',ylim=range(c(f.pred[,i])),xlab=xlab.name,ylab=ylab.name)
      }
      if (U.ps[i]=='nonhomo'){
        order.temp=order(U.pred[,i])
        plot(U.pred[,i][order.temp],f.pred[,i][order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
      }
    }
  }
    return (list(beta=beta.hat,f=f.hat,f.pred=f.pred,record.beta=record.beta,record.h=record.h))
}


fit.aplm.bootstrap<-function(Z,U,X,original.Weight,PoolID,U.ps,X.pred,U.pred,kernel_type,max.iter=10,plot=FALSE,plot.pred=FALSE,BS=50){
  
  true.Z=as.matrix(Z)
  true.U=as.matrix(U)
  true.X=as.matrix(X)
  true.PoolID=PoolID
  true.original.Weight=original.Weight
  
  bs.res=list()
  
  for (bs in 1:BS){
    
    
    bs.ID.list=sample(unique(true.PoolID),replace = TRUE)
    
    
    bs.Z=matrix(nrow=0,ncol=ncol(true.Z))
    bs.U=matrix(nrow=0,ncol=ncol(true.U))
    bs.X=matrix(nrow=0,ncol=ncol(true.X))
    bs.original.Weight=c()
    bs.PoolID=c()
    bs.count=0
    for (bs.ID in bs.ID.list){
      bs.count=bs.count+1
      index=which(true.PoolID==bs.ID)
      bs.Z=rbind(bs.Z,as.matrix(true.Z[index,]))
      bs.X=rbind(bs.X,as.matrix(true.X[index,]))
      bs.U=rbind(bs.U,as.matrix(true.U[index,]))
      bs.original.Weight=c(bs.original.Weight,true.original.Weight[index])
      bs.PoolID=c(bs.PoolID,rep(bs.count,length(index)))
    }
    
    Z=bs.Z
    X=bs.X
    U=bs.U
    PoolID=bs.PoolID
    original.Weight=bs.original.Weight  
    
    bs.data=as.data.frame(cbind(Z,U,X,PoolID,original.Weight))
    
    
    Z=as.matrix(Z)
    U=as.matrix(U)
    X=as.matrix(X)
    
    p=ncol(X)
    q=ncol(U)
    N=nrow(Z)
    J=length(unique(PoolID))
    
    #sort all data by PoolID
    index=order(PoolID)
    Z=Z[index]
    U=U[index,]
    X=X[index,]
    original.Weight=original.Weight[index]
    PoolID=PoolID[index]
    
    #identify the start and end location of each pool
    unique_PoolID=unique(PoolID)
    start_ind=rep(0,length(unique_PoolID))
    end_ind=rep(0,length(unique_PoolID))
    count=0
    for (id in unique_PoolID){
      count=count+1
      index=which(PoolID==id)
      start_ind[count]=index[1]
      end_ind[count]=index[length(index)]
    }
    
    Weight=rep(NA,N)
    for (i in 1:length(start_ind)){
      index=start_ind[i]:end_ind[i]
      Weight[index]=original.Weight[index]/sum(original.Weight[index])
    }
    
    #pre calculate some variables to speed up
    #pre calculation for linear
    # H=solve(t(X)%*%X)%*%t(X)
    X.pool=matrix(NA,nrow=J,ncol=p)
    Z.pool=rep(NA,J)
    original.Weight.pool=rep(NA,J)
    for (j in 1:J){
      index=start_ind[j]:end_ind[j]
      X.pool[j,]=colSums(matrix(X[index,]*Weight[index],nrow=length(index),ncol=p))
      Z.pool[j]=Z[index][1]
      original.Weight.pool[j]=sum(original.Weight[index])
    }
    H.pool=solve(t(X.pool*original.Weight.pool)%*%X.pool)%*%t(X.pool*original.Weight.pool)
    
    #pre calculation for nonlinear
    U.pool=matrix(NA,ncol=q,nrow=J)
    for (j in 1:J){
      index=start_ind[j]:end_ind[j]
      for (i in 1:q){
        U.pool[j,i]=sum(U[index,i]*Weight[index])
      }
    }
    
    # #pre define vectors to record results
    record.beta=matrix(NA,nrow=p,ncol=max.iter)
    record.f=array(0,c(N,q,max.iter))
    record.h=matrix(NA,nrow=q,ncol=max.iter)
    
    #initialize
    mu.hat=mean(Z)
    f.hat=matrix(0,nrow=N,ncol=q)
    Y.hat.pre=rep(mu.hat,N)
    beta.hat.pre=rep(0,p)
    cond=TRUE
    iter=0
    
    #algorithm starts
    while (iter==0 | cond==TRUE){
      iter=iter+1
      Y.hat=matrix(rep(NA,N),nrow=N,ncol=1)
      for (j in 1:J){
        index=start_ind[j]:end_ind[j]
        for (i in 1:length(index)){
          Y.hat[index[i]]=(Z[index[i]]-sum(Weight[index[-i]]*Y.hat.pre[index[-i]]))/Weight[index[i]]
        }
      }
      f.hat.pool=matrix(NA,nrow=J,ncol=q)
      for (j in 1:J){
        index=start_ind[j]:end_ind[j]
        f.hat.pool[j,]=colSums(matrix(f.hat[index,]*Weight[index],nrow=length(index),ncol=q))
      } 
      beta.hat=H.pool%*%(Z.pool-rowSums(f.hat.pool))
      # beta.hat=coef.beta
      X.beta.hat=X%*%beta.hat

      #calculate X.beta.hat.pool in preparation for homogenous pooling in nonlinear part
      X.beta.hat.pool=rep(NA,J)
      for (j in 1:J){
        index=start_ind[j]:end_ind[j]
        X.beta.hat.pool[j]=sum(X.beta.hat[index]*Weight[index])  
      }
      
      #update nonlinear part
      for (i in 1:q){
        if (U.ps[i]=='nonhomo'){ #random pooling Ver.
          cv_h<-function(h){
            Y.temp=Y.hat-X.beta.hat-rowSums(as.matrix(f.hat[,-i]))
            res=rep(NA,N)
            for (j in 1:J){
              index=start_ind[j]:end_ind[j]
              Sc.temp=LLS_W_cpp(U[-index,i],h,U[index,i],original.Weight[-index],kernel_type) #S1
              res[index]=(Y.temp[index]-Sc.temp%*%Y.temp[-index])^2
            }
            return (mean(res))
          }
          tt=proc.time()
          print(paste("The ",iter,'th iteration, optimizing h',i,'.',sep=''))
          h=optimize(cv_h,lower=diff(range(U[,i]))/5,upper=diff(range(U[,i]))/2)$minimum
          print(paste('Done. h',i,'=',h,' Time elapse is ',round((proc.time()-tt)[3],2),sep=''))
          S1=LLS_W_cpp(U[,i],h,U[,i],original.Weight,kernel_type)
          S1c=S1-matrix(colMeans(S1),byrow=TRUE,nrow=N,ncol=N)
          f.hat[,i]=S1c%*%(Y.hat-X.beta.hat-rowSums(as.matrix(f.hat[,-i])))
        }else{  #homopooling Ver.

          f.hat.pool=matrix(NA,nrow=J,ncol=q)
          for (j in 1:J){
            index=start_ind[j]:end_ind[j]
            f.hat.pool[j,-i]=colSums(matrix(f.hat[index,-i]*Weight[index],nrow=length(index),ncol=q-1))
          }
          cv_h<-function(h){
            Y.temp=Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool
            res=rep(NA,J)
            for (j in 1:J){
              Sc.temp=LLS_W_cpp(U.pool[-j,i],h,U.pool[j,i],original.Weight.pool[-j],kernel_type) 
              res[j]=(Y.temp[j]-Sc.temp%*%Y.temp[-j])^2
            }
            return (mean(res))
          }
          tt=proc.time()
          print(paste("The ",iter,'th iteration, optimizing h',i,'.',sep=''))
          h=optimize(cv_h,lower=diff(range(U[,i]))/5,upper=diff(range(U[,i]))/2)$minimum
          print(paste('Done. h',i,'=',h,' Time elapse is ',round((proc.time()-tt)[3],2),sep=''))
          S1=LLS_W_cpp(U.pool[,i],h,U.pool[,i],original.Weight.pool,kernel_type)
          S1c=S1-matrix(colMeans(S1),byrow=TRUE,nrow=J,ncol=J)
          f.hat.pool[,i]=S1c%*%(Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool)  
          for (j in 1:J){
            index=start_ind[j]:end_ind[j]
            f.hat[index,i]=f.hat.pool[j,i]
          }
        }
        record.h[i,iter]=h
      }
      
      Y.hat.pre=X.beta.hat+rowSums(f.hat)
      cond= (iter<max.iter)
      beta.hat.pre=beta.hat
      
      
      #record the result of each iteration
      record.beta[,iter]=beta.hat
      record.f[,,iter]=f.hat
      # record.h has been saved
      
      #plot part
      if (plot==TRUE){
        par(mfrow=c(1,q))
        sqerror=0
        for (i in 1:q){
          xlab.name=expression(('z')[i])
          ylab.name=expression(hat('f')[i]*'(z)')
          if (U.ps[i]=='homo'){
            Y.plot.temp=rep(NA,J)
            U.plot.temp=rep(NA,J)
            for (j in 1:J){
              index=start_ind[j]:end_ind[j]
              Y.plot.temp[j]=sum(Weight[index]*f.hat[index,i])
              U.plot.temp[j]=sum(Weight[index]*U[index,i])
            }
            order.temp=order(U.plot.temp)
            plot(U.plot.temp[order.temp],Y.plot.temp[order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
          }else{
            order.temp=order(U[,i])
            plot(U[,i][order.temp],f.hat[,i][order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
          }
        }
      }
      Sys.sleep(1)
    }
    #delete unused results due to earlier stop than max.iter (reached accuracy requirement)
    record.beta=as.matrix(record.beta[,1:iter])
    record.f=as.matrix(record.f[,,1:iter])
    record.h=as.matrix(record.h[,1:iter])
    
    #prediction part
    f.pred=matrix(NA,nrow=nrow(as.matrix(U.pred)),ncol=q)
    for (i in 1:q){
      if (U.ps[i]=='nonhomo'){
        S1.pred=LLS_W_cpp(U[,i],record.h[i,iter],U.pred[,i],original.Weight,kernel_type)
        S1.original=LLS_W_cpp(U[,i],record.h[i,iter],U[,i],original.Weight,kernel_type)
        S1c=S1.pred-matrix(colMeans(S1.original),byrow=TRUE,nrow=nrow(f.pred),ncol=N)
        f.pred[,i]=S1c%*%(Y.hat.pre-X.beta.hat-rowSums(as.matrix(f.hat[,-i]))) #use Y.hat.pre
      }
      if (U.ps[i]=='homo'){
        #predict based on homo pooled results
        S1.pred=LLS_W_cpp(U.pool[,i],record.h[i,iter],U.pred[,i],original.Weight.pool,kernel_type)
        S1.original=LLS_W_cpp(U.pool[,i],record.h[i,iter],U.pool[,i],original.Weight.pool,kernel_type)
        S1c=S1.pred-matrix(colMeans(S1.original),byrow=TRUE,nrow=nrow(f.pred),ncol=J)
        f.pred[,i]=S1c%*%(Z.pool-rowSums(as.matrix(f.hat.pool[,-i]))-X.beta.hat.pool)
      }
    }
    #Y.pred
    Y.pred=X.pred%*%beta.hat+rowSums(as.matrix(f.pred))
    if (plot.pred==TRUE){
      par(mfrow=c(1,q))
      sqerror=0
      for (i in 1:q){
        xlab.name=paste('z',i,sep='')
        ylab.name=paste('f',i,'(z',i,')',sep='')
        if (U.ps[i]=='homo'){
          Y.plot.temp=(f.pred[,i])
          U.plot.temp=(U.pred[,i])
          
          order.temp=order(U.plot.temp)
          plot(U.plot.temp[order.temp],Y.plot.temp[order.temp],type='l',ylim=range(c(f.pred[,i])),xlab=xlab.name,ylab=ylab.name)
        }
        if (U.ps[i]=='nonhomo'){
          order.temp=order(U.pred[,i])
          plot(U.pred[,i][order.temp],f.pred[,i][order.temp],type='l',xlab=xlab.name,ylab=ylab.name)
        }
      }
    }
    bs.res[[bs]]=list(beta=beta.hat,f=f.hat,f.pred=f.pred,record.beta=record.beta,record.h=record.h)
  }
  
  return (bs.res)
}

