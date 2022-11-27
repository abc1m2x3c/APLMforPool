#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Gaussian kernel 
// [[Rcpp::export]]
vec GKernel(vec x){
  vec res(x.size());
  res=exp(-pow(x,2)/2);
  return(res);
}

// Ep kernel 
// [[Rcpp::export]]
vec EKernel(vec x){
  vec res(x.size());
  res=(1-pow(x,2))*3/4;
  res.elem(find(pow(x,2)>1)).zeros();
  return(res);
}


//return a weight vector for each nx
// [[Rcpp::export]]
mat LLS_W_cpp(vec U, double h, vec nx , vec Weight,int kernel){
  vec Uu(U.size());
  vec W_u(U.size());
  mat res(nx.size(),U.size());
  mat Gamma_u(U.size(),2,fill::ones);
  mat W_u_Gamma_u(size(Gamma_u));
  for(unsigned int k=0; k<nx.size(); k++){
    double u=nx[k];
    Uu=U-u;
    Gamma_u.col(1)=Uu;
    if(kernel==1) //Gaussian kernel
    {
      W_u=1/h*GKernel(Uu/h)%Weight;
    }
    if(kernel==2) //Ep kernel
    {
      W_u=1/h*EKernel(Uu/h)%Weight;
    }
    // Gamma_u=join_cols(X,Uu_X);
    W_u_Gamma_u=Gamma_u;
    W_u_Gamma_u.each_col() %= W_u;
    mat temp(2,U.size(),fill::zeros);
    //temp=((Gamma_u.t()*W_u_Gamma_u).i()*(W_u_Gamma_u.t())); //only differece with LL_cpp
    temp=solve((Gamma_u.t()*W_u_Gamma_u),(W_u_Gamma_u.t())); //only differece with LL_cpp
    res.row(k) = temp.row(0);
  }
  return (res);
}


//return eta(nx) for each nx; have some errors due to Y
// [[Rcpp::export]]
mat LL_W_cpp(vec U, vec Y,double h, vec nx , vec Weight,int kernel){
  vec Uu(U.size());
  vec W_u(U.size());
  mat res(nx.size(),U.size());
  mat Gamma_u(U.size(),2,fill::ones);
  mat W_u_Gamma_u(size(Gamma_u));
  for(unsigned int k=0; k<nx.size(); k++){
    double u=nx[k];
    Uu=U-u;
    Gamma_u.col(1)=Uu;
    if(kernel==1) //Gaussian kernel
    {
      W_u=1/h*GKernel(Uu/h)%Weight;
    }
    if(kernel==2) //Ep kernel
    {
      W_u=1/h*EKernel(Uu/h)%Weight;
    }
    // Gamma_u=join_cols(X,Uu_X);
    W_u_Gamma_u=Gamma_u;
    W_u_Gamma_u.each_col() %= W_u;
    mat temp(2,U.size(),fill::zeros);
    temp=(Gamma_u.t()*W_u_Gamma_u).i()*W_u_Gamma_u.t(); //only differece with LL_cpp
    res.row(k) = temp.row(0);
  }
  return (res);
}


// [[Rcpp::export]]
mat deleRows(mat x, int k1, int k2){
  x.shed_rows(k1-1,k2-1);
  return(x);
}


//; have some errors due to LL_W_cpp
// [[Rcpp::export]]
double LL_W_cv_cpp(vec U, vec Y, vec Weight, double h, IntegerVector start_ind, IntegerVector end_ind, int kernel){
  vec res(start_ind.size());
  for (int index=0; index<start_ind.size(); index++){
    vec u_grid(end_ind(index)-start_ind(index)+1);
    u_grid=U.rows(start_ind[index]-1,end_ind[index]-1);
    // mat X_rest=deleRows(X,start_ind[index],end_ind[index]);
    vec U_rest=deleRows(U,start_ind[index],end_ind[index]);
    vec Y_rest=deleRows(Y,start_ind[index],end_ind[index]);
    vec Weight_rest=deleRows(Weight,start_ind[index],end_ind[index]);
    vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
    res[index]=sum(pow(Y_keep-LL_W_cpp(U_rest,Y_rest, h,u_grid , Weight_rest,kernel),2));
    
    // mat beta = VCM_cpp(X_rest, U_rest, Y_rest, h, u_grid, kernel);
    // beta.shed_rows(p,2*p-1);
    // mat X_keep=X.rows(start_ind[index]-1,end_ind[index]-1);
    // vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
    // res[index]=sum(pow(Y_keep-sum(X_keep%beta.t(),1),2));
  }
  return (mean(res));
}


// // [[Rcpp::export]]
// double VCM_cv_cpp(mat X, vec U, vec Y, double h, vec W_sample, vec W_cv,IntegerVector start_ind, IntegerVector end_ind, int kernel){
//   int p = X.n_cols;
//   vec res(start_ind.size());
//   for (int index=0; index<start_ind.size(); index++){
//     vec u_grid(end_ind(index)-start_ind(index)+1);
//     u_grid=U.rows(start_ind[index]-1,end_ind[index]-1);
//     mat X_rest=deleRows(X,start_ind[index],end_ind[index]);
//     vec U_rest=deleRows(U,start_ind[index],end_ind[index]);
//     vec Y_rest=deleRows(Y,start_ind[index],end_ind[index]);
//     vec W_sample_rest=deleRows(W_sample,start_ind[index],end_ind[index]);
//     mat beta = VCM_cpp(X_rest, U_rest, Y_rest, h, u_grid, W_sample_rest, kernel);
//     beta.shed_rows(p,2*p-1);
//     mat X_keep=X.rows(start_ind[index]-1,end_ind[index]-1);
//     vec Y_keep=Y.rows(start_ind[index]-1,end_ind[index]-1);
//     vec W_cv_keep=W_cv.rows(start_ind[index]-1,end_ind[index]-1);
//     res[index]=sum(W_cv_keep%pow(Y_keep-sum(X_keep%beta.t(),1),2));
//   }
//   return (mean(res));
// }


