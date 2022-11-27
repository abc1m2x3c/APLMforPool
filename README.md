# APLMforPool
This project is to fit APLM for pooled biomonitoring data. Paper is currently under review at Computational Statistics & Data Analysis.

The necessary source code is src.R and PLM.cpp. Please download it to your working directory of R or Rstudio.

## A simple example 
Below is an simple example to use the codes. All codes below are included in 'Example.R' in this repository. It analyzed the PBDE-47 concentration collected by NHANES during 2015-2016. You can download the dataset from NHANES official website. We have prepared a cleaned dataset. You can find it in 'PBDE47-2015-2016.csv' in this repository.


### Step 1. Clean the memory, install (if has not) and load the required R packages, and source the code in this repository  
'''
packages <- c("maxLik","Rcpp","RcppArmadillo")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
library(Rcpp)
library(RcppArmadillo)
source('src.R')
'''
