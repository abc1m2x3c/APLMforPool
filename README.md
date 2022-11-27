# APLMforPool
This project is to fit APLM for pooled biomonitoring data. Paper is currently under review at Computational Statistics & Data Analysis.

The necessary source code is src.R and PLM.cpp. Please download it to your working directory of R or Rstudio.

## A simple example 
Below is an simple example to use the codes. All codes below are included in 'Example.R' in this repository. It analyzed the PBDE-47 concentration collected by NHANES during 2015-2016. You can download the dataset from NHANES official website. We have prepared a cleaned dataset. You can find it in 'PBDE47-2015-2016.csv' in this repository.

### Step 1. Clean the memory, install (if has not) and load the required R packages, and source the code in this repository  
```
packages <- c("maxLik","Rcpp","RcppArmadillo")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
library(Rcpp)
library(RcppArmadillo)
source('src.R')
```

### Step 2. Input data
```
input_data=read.csv('PBDE47-2015-2016.csv')
```
The dataset should be in the following form, where the first column Z is the pooled concentration levels
![Optional Text](https://github.com/abc1m2x3c/APLMforPool/blob/f9a26ba36a96c11df126a8c39a8221cefd9094e6/DatasetScreenshot.PNG)

### Step 3. Determine the input argument

```
kernel_type=2 #kernel_type=1 Gaussian kernel; 2 Epanechnikov kernel
Z=input_data$Z # pooled concentration
U=cbind(log(input_data$BMI),input_data$age) #covariates in nonlinear component
                                            #In this example, U1 is log(NMI), U2 is age
X=cbind(1,input_data$race,input_data$gender) #covariates in linear component; the first column is a vector of 1's for intercept
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
```

### Step 4. Fit the model
```
res=fit.aplm(Z,U,X,original.Weight,PoolID,U.ps,X.pred=X.pred,U.pred=U.pred,kernel_type=kernel_type,max.iter=max.iter,plot=TRUE,plot.pred=TRUE)
```
The program will output the selected bandwidth in each iteration. If plot=TRUE and plot.pred=TRUE are specified, it will also output the fitted curves in each iteration and the predicted curves in the end of the program, as shown below:


