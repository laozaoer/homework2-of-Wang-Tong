#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat SVT(arma::mat X,double lambda){
  arma::mat U,V,Ud,Vd,Xhat;arma::vec s,sd;
  arma::svd(U,s,V,X);
  Ud=U.cols(find(s>0.000001));
  Vd=V.cols(find(s>0.000001));
  sd=s.elem(find(s>0.000001));
  for(int i=0;i<sd.n_elem;i++){
    if(sd(i)>lambda){sd(i)=sd(i)-lambda;}
    else{sd(i)=0;}
  }
  Xhat=Ud*diagmat(sd)*trans(Vd);
  return Xhat;
}
arma::mat shedmat0(arma::mat a,double p){
  int m;int n;
  m=a.n_rows;n=a.n_cols;
  arma::mat b;b.randu(m,n);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      if(b(i,j)<=p){a(i,j)=0;}
    }
  }
  return a;
}
arma::mat projectmat(arma::mat a,arma::mat b){
  for(int i=0;i<b.n_rows;i++){
    for(int j=0;j<b.n_cols;j++){
      if(b(i,j)==0){a(i,j)=0;}
    }
  }
  return a;
}
// [[Rcpp::export]]
arma::mat shedmat1(arma::mat a,double p){
  int m;int n;
  m=a.n_rows;n=a.n_cols;
  arma::mat b;b.randu(m,n);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      if(b(i,j)<=p){a(i,j)=1;}
    }
  }
  return a;
}
// [[Rcpp::export]]
arma::mat APO2(arma::mat X,double lambda,double p,double c,double beta){
  arma::mat Xhat,Xmiss,Yhat,Zhat,tempX;int m,n;m=X.n_rows;n=X.n_cols;double delta=1;
  Xmiss=shedmat0(X,p);Zhat.randu(m,n);
  Xhat.randu(m,n);
  do{
    Yhat=Zhat+delta*projectmat((Xmiss-Zhat),Xmiss);
    tempX=Xhat;
    Xhat=SVT(Yhat,lambda*delta);
    Zhat=Xhat+beta*(Xhat-tempX);
  } while (norm(tempX-Xhat)>c);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      if(Xhat(i,j)>1){Xhat(i,j)=1;}
      if(Xhat(i,j)<0){Xhat(i,j)=0;}
    }
  }
  return Xhat;
}