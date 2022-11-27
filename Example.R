packages <- c("maxLik","Rcpp","RcppArmadillo")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
library(Rcpp)
library(RcppArmadillo)
source('src.R')

input_data=read.csv('PBDE47-2015-2016.csv')
kernel_type=2 #kernel_type=1 Gaussian kernel; 2 Epanechnikov kernel
Z=input_data$Z # pooled concentration
U=cbind(log(input_data$BMI),input_data$age) #covariates in nonlinear component
                                            #In this example, U1 is log(NMI), U2 is age
X=cbind(1,input_data$race,input_data$gender) #covariates in linear component
                                             #In this example, X1 is race with 1 white 2 non-white
                                             #                 X2 is gender with 1 male 2 female   
original.Weight=input_data$Weight #sampling weight of each individual
PoolID=input_data$PoolID #Pool ID
U.ps=c('nonhomo','homo') #pooling structure corresponding to U
                         #In this example, BMI is non-homogeneous pooling and age is homogeneous pooling
max.iter=10 #maximum iteration
len=100 #length of covariates for prediction
U.pred=cbind((seq(min(U[,1]),max(U[,1]),length=len)),(seq(min(U[,2]),max(U[,2]),length=len))) # U for prediction
X.pred=matrix(c(1,rep(0,ncol(X)-1)),nrow=len,ncol=ncol(X),byrow = TRUE) # X for prediction

#fit the model
#output includes:
# beta:   final estimation of beta
# f:      an N by q matrix with column 1 to q the estimation of f1 to fq at each U
# f.pred: a len by q matrix with column 1 to q the estimation of f1 to fq at the corresponding U.pred
# record.h: recording the bandwidths at each iteration
res=fit.aplm(Z,U,X,original.Weight,PoolID,U.ps,X.pred=X.pred,U.pred=U.pred,kernel_type=kernel_type,max.iter=max.iter,plot=TRUE,plot.pred=TRUE)

res$beta
plot(U.pred[,1],res$f.pred[1:len,1],col="black",ylab='',type='l',xlab='log(BMI)',cex.lab=1.5,cex.axis=1.1,lwd=2)
title(ylab=expression(hat('f')[1]*'(log(BMI))'), line=2, cex.lab=1.5)
plot(U.pred[,2],res$f.pred[1:len,2],col="black",ylab='',type='l',xlab='Age',cex.lab=1.5,cex.axis=1.1,lwd=2)
title(ylab=expression(hat('f')[2]*'(age)'), line=2, cex.lab=1.5)


