packages <- c("maxLik","Rcpp","RcppArmadillo")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
library(Rcpp)
library(RcppArmadillo)
source('src.R')

input_data=read.csv('PBDE47-2015-2016.csv')

kernel_type=2
Z=input_data$Z
U=cbind(log(input_data$BMI),input_data$age)
X=cbind(1,input_data$race,input_data$gender)
original.Weight=input_data$Weight
PoolID=input_data$PoolID
U.ps=c('nonhomo','homo')
max.iter=10
len=100
U.pred=cbind((seq(min(U[,1]),max(U[,1]),length=len)),(seq(min(U[,2]),max(U[,2]),length=len)))
X.pred=matrix(c(1,rep(0,ncol(X)-1)),nrow=len,ncol=ncol(X),byrow = TRUE)

res=fit.aplm(Z,U,X,original.Weight,PoolID,U.ps,X.pred=X.pred,U.pred=U.pred,kernel_type=kernel_type,max.iter=max.iter,plot=TRUE,plot.pred=TRUE)


plot(U.pred[,1],res$f.pred[1:len,1],col="black",ylab='',type='l',xlab='log(BMI)',cex.lab=1.5,cex.axis=1.1,lwd=2)
title(ylab=expression(hat('f')[1]*'(log(BMI))'), line=2, cex.lab=1.5)
plot(U.pred[,2],res$f.pred[1:len,2],col="black",ylab='',type='l',xlab='Age',cex.lab=1.5,cex.axis=1.1,lwd=2)
title(ylab=expression(hat('f')[2]*'(age)'), line=2, cex.lab=1.5)

# 
# 
# BS=2
# res.bs=fit.aplm.bootstrap(Z,U,X,original.Weight,PoolID,U.ps,X.pred=X.pred,U.pred=U.pred,kernel_type=kernel_type,max.iter=max.iter,plot=FALSE,plot.pred=FALSE,BS=BS)
# filename=paste(chem,'maxiter',max.iter,'BS',BS,'ID',SLURM_ID,'.RData',sep='')
# save(res.bs,file=filename)
# 
