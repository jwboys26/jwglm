#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>



const double pi = arma::datum::pi;
const double M_sqrt2pi=2.506628;



double normalCDF(double value){
  double out=0;
  const double M_sqrt2=1.414214;
  
  out = 0.5 * erfc(-value / M_sqrt2);
  return out;
}



double fl(double x, Rcpp::String strDistr="Logit"){
  
  double out = 0;
  double ex = exp(x);
  
  if(strDistr=="Probit"){
    out = exp(-pow(x,2)/2)/M_sqrt2pi;
  }else if(strDistr=="LogLog"){
    out = ex*exp(-ex);
  }else{
    
    out = exp(-x)/pow( (1+exp(-x)), 2 );
  }
  
  return out;
  
}


double fl_dot(double x, Rcpp::String strDistr="Logit"){
  
  double out = 0;
  double ex = exp(x);
  double emx = exp(-x);
  
  
  if(strDistr=="Probit"){
    out = -x*exp(-pow(x,2)/2)/M_sqrt2pi;
  }else if(strDistr=="LogLog"){
    out = ex*exp(-ex)*(1-ex);
  }else{
    
    out = -emx *(1-emx)/pow( (1+emx),3 ) ;
  }
  
  return out;
  
}

double fl_2dot(double x, Rcpp::String strDistr="Logit"){
  
  double out = 0;
  double ex = exp(x);
  double emx = exp(-x);
  double emx2 = exp(-2*x);
  
  if(strDistr=="Probit"){
    out = (pow(x,2)-1)*exp(-pow(x,2)/2)/M_sqrt2pi;
  }else if(strDistr=="LogLog"){
    out = ex*exp(-ex)*(1-3*ex+exp(2*x));
  }else{
    
    out = -emx2 *(4-ex-emx)/pow( (1+emx),4 ) ;
  }
  
  return out;
  
}


//[[Rcpp::export]]
double Fl(double x, Rcpp::String strDistr="Logit"){
  
  double out = 0;
  double ex = exp(x);
  
  if(strDistr =="Probit"){
    out = normalCDF(x);
  }else if(strDistr=="LogLog"){
    out = 1-exp(-ex);
  }else{
    
    out = 1/ ( 1+exp(-x) );
  }
  
  return out;
  
}



double GetVdj(arma::vec Y, arma::mat D, int j){
  
  int n = Y.n_elem;
  double ans = 0;
  
  int indval=0;
  int Yk=0;
  double dkj=0;
  
  for(int k=1; k<=n; k++){
    
    
    Yk = Y[k-1];
    dkj = D(k-1, j-1);
    
    if(Yk==1){
      indval = 1;
    }else{
      indval = 0;
    }
    
    ans += dkj*indval;
    
  }
  
  return ans;
}



double GetJdj(arma::mat D, int j, arma::vec pVec, Rcpp::String strDistr){
  
  int n = pVec.n_elem;
  double out = 0;
  double pk = 0;
  double dkj=0;
  
  
  for(int k=1;k<=n;k++){
    pk = pVec[k-1];
    dkj = D(k-1, j-1);
    out += dkj*pk; 
  }  
    
  return out;
}



double GetWdj(arma::vec Y, arma::mat D, int j, arma::vec pVec, Rcpp::String strDistr){
  
  double out=GetVdj(Y, D, j) - GetJdj(D, j, pVec, strDistr);
  return out;
}

arma::vec Get_pVec(arma::mat X, arma::vec beta, Rcpp::String strDistr){
  
  int n = X.n_rows;
  arma::vec out(n);
  out.zeros();
  
  arma::vec Xb = X*beta;
  double x = 0;
  
  for(int i=1;i<=n;i++){
    x = Xb[i-1];
    out[i-1] = Fl(x, strDistr);
  }
  
  return out;
  
}




arma::vec GetSj(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, int j, arma::vec pVec, Rcpp::String strDistr){
  
  int n = Y.n_elem;
  int p = X.n_cols;
  
  double Wdj=0;
  double dkj = 0;
  
  
  arma::rowvec xk(p);
  arma::vec xkb(1);
  double pkb = 0;
  
  arma::vec out(p);
  out.zeros();
  
  
  Wdj = GetWdj(Y, D, j, pVec, strDistr);
  for(int k=1;k<=n;k++){
    xk = X.row(k-1);
    
    xkb = xk * beta;
    pkb = fl(xkb[0]);
    
    dkj = D(k-1,j-1);
    out += Wdj* dkj*pkb * xk.t();
    
  }
  
  
  out *= -2;
  
  return out;
  
}

arma::vec GetSn(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, arma::vec pVec, Rcpp::String strDistr){
  
  int p = X.n_cols;
  arma::vec out(p);
  out.zeros();
  
  for(int j=1;j<=p; j++){
    
    out += GetSj(beta, Y, X, D, j, pVec, strDistr);
  }
  
  return out;
  
}


arma::mat GetWnj(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, int j, arma::vec pVec, Rcpp::String strDistr){
  
  int n = Y.n_elem;
  int p = X.n_cols;
  

  double dkj = 0;
  double dij = 0;
  
  
  arma::rowvec xk(p);
  arma::vec xkb(1);
  double pkb = 0;
  
  arma::rowvec xi(p);
  arma::vec xib(1);
  double pib = 0;
  
  
  arma::mat out(p, p);
  out.zeros();
  
  
  for(int i=1;i<=n;i++){
    
    xi = X.row(i-1);
    
    xib = xi * beta;
    pib = fl(xib[0]);
    
    dij = D(i-1,j-1);
    
    
    for(int k=1;k<=n;k++){
      xk = X.row(k-1);
      
      xkb = xk * beta;
      pkb = fl(xkb[0]);
      
      dkj = D(k-1,j-1);
      out += dij*pib * dkj*pkb * xk.t()*xi;
      
      
    }
    
    
  }
  
  return out;
  
  
}


arma::mat GetWn(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, arma::vec pVec, Rcpp::String strDistr){
  
  int p = X.n_cols;
  arma::mat out(p,p);
  out.zeros();
  
  for(int j=1;j<=p; j++){
    
    out += GetWnj(beta, Y, X, D, j, pVec, strDistr);
  }
  
  return out;
  
  
}
  



double GetLp(arma::vec Y, arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr){
  
  arma::vec pVec = Get_pVec(X, beta, strDistr);
  
  double out=0;
  
  int p = D.n_cols;
  
  double val = 0;
  
  for(int j=1; j<=p;j++){
    val = GetWdj(Y,D, j, pVec, strDistr);
    out += pow(val,2);
    
  }
  
  return out;
}


// [[Rcpp::export]]
arma::vec Get_beta_only(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, 
              Rcpp::String strDistr, int nIter=100, 
              double lr=0.01, double crit=1e-3, bool bDisp=false){
  
  int p = X.n_cols;
  arma::vec grad(p);
  grad.zeros();
  
  arma::vec newbeta(p);
  newbeta.zeros();
  
  arma::vec oldbeta(p);
  oldbeta.zeros();
  
  arma::vec diffvec(p);
  diffvec.zeros();
  
  
  oldbeta = beta;
  
  double diff = 0;
  
  arma::vec pVec(p);
  
  int nFlag = 0;
  double lossval=0;
  
  arma::vec LossVec(nIter);
  
  arma::vec Sn(p);
  Sn.zeros();
  
  for(int iter=1;iter<=nIter;iter++){
    
    pVec = Get_pVec(X, oldbeta, strDistr);
    Sn = GetSn(oldbeta, Y, X, D, pVec, strDistr);
    
    newbeta = oldbeta - lr*Sn; 
    
    
    diffvec = newbeta-oldbeta;
    diff = sqrt(sum(diffvec%diffvec));
    
    lossval = GetLp(Y, X, D, newbeta, strDistr);
    LossVec[iter-1] = lossval;
    
    
    if(diff<crit){
      
      LossVec = LossVec.subvec(0, iter-1);
      nFlag = iter;
      
      break;
      
    }
    oldbeta = newbeta;
    nFlag = iter;
    
  }
  
  if( (nFlag==nIter) & (bDisp==true)){
    
    Rcout<< "The convergence of beta was not reached." << std::endl;
  }
  
  arma::vec out = newbeta;

  return out;
  
  
}




List Get_beta(arma::vec beta, arma::vec Y, arma::mat X, arma::mat D, 
                Rcpp::String strDistr, int nIter=100, 
                double lr=0.01, double crit=1e-3, bool bDisp=false){
  
  int p = X.n_cols;
  arma::vec grad(p);
  grad.zeros();
  
  arma::vec newbeta(p);
  newbeta.zeros();
  
  arma::vec oldbeta(p);
  oldbeta.zeros();
  
  arma::vec diffvec(p);
  diffvec.zeros();
  
  
  oldbeta = beta;
  
  double diff = 0;
  
  arma::vec pVec(p);
  
  int nFlag = 0;
  double lossval=0;
  
  arma::vec LossVec(nIter);
  
  arma::vec Sn(p);
  Sn.zeros();
  
  for(int iter=1;iter<=nIter;iter++){
    
    pVec = Get_pVec(X, oldbeta, strDistr);
    Sn = GetSn(oldbeta, Y, X, D, pVec, strDistr);
    
    newbeta = oldbeta - lr*Sn; 
    

    diffvec = newbeta-oldbeta;
    diff = sqrt(sum(diffvec%diffvec));
    
    lossval = GetLp(Y, X, D, newbeta, strDistr);
    LossVec[iter-1] = lossval;
    
    
    if(diff<crit){
      
      LossVec = LossVec.subvec(0, iter-1);
      nFlag = iter;
      
      break;
      
    }
    oldbeta = newbeta;
    nFlag = iter;
    
  }
  
  if( (nFlag==nIter) & (bDisp==true)){
  
    Rcout<< "The convergence of beta was not reached." << std::endl;
  }
  
  
  List lst(4);
  lst[0] = nFlag;
  lst[1] = diff;
  lst[2] = newbeta;
  lst[3] = LossVec;
  
  

  return lst;
  
  
}


arma::mat Get_Wn(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr="Logit"){
  
  int n = X.n_rows;
  int J = X.n_cols;
  
  
  arma::mat out(J,J);
  out.zeros();
  double dij, dkj;
  double fi, fk;
  
  arma::mat tmp(1,1);
  arma::rowvec qi(J);
  qi.zeros();
  
  arma::rowvec qk(J);
  qk.zeros();
  
  
  
  for(int j=1;j<=J;j++){
    
    for(int i=1;i<=n;i++){
      
      dij = D(i-1,j-1);
      tmp = X.row(i-1)*beta;
      
      fi = fl(tmp(0,0), strDistr=strDistr);
      qi = fi*X.row(i-1);
      
      for(int k=1;k<=n;k++){
        dkj = D(k-1,j-1);
        tmp = X.row(k-1)*beta;
        
        fk = fl(tmp(0,0), strDistr=strDistr);
        qk = fk*X.row(k-1);
        
        out += dij*dkj*qk.t()*qi;
      }
      
      
    }
    
  }
  
  
  out *=2;
  
  return out;
  
  
}





arma::vec Get_d(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr="Logit"){
  
  int n = X.n_rows;
  int J = X.n_cols;
  
  arma::mat Wn = Get_Wn(X, D, beta);
  
  arma::mat Pn(n,n);
  Pn.zeros();
  
  arma::mat Lambda_n(n,n);
  Lambda_n.zeros();
  
  arma::mat tmp(1,1);
  double Fk=0;
  
  for(int k=1;k<=n;k++){
    tmp = X.row(k-1)*beta;
    Fk = Fl(tmp(0,0), strDistr=strDistr);
    Pn(k-1,k-1) = Fk*(1-Fk);
    Lambda_n(k-1,k-1) = fl(tmp(0,0), strDistr=strDistr);
  }
  
  arma::mat P2 = D*D.t()*Pn*D*D.t();
  arma::mat Wtilde = inv(Wn)*X.t()*Lambda_n;
  
  arma::mat Wtilde2(n*J, pow(n,2));
  Wtilde2.zeros();
  
  int sidx=0, eidx=0;
  int sidx2=0, eidx2=0;
  
  for(int i=1;i<=n;i++){
    sidx = (i-1)*J+1;
    eidx = sidx+J-1;
    
    sidx2 = (i-1)*n+1;
    eidx2 = sidx2+n-1;
    
    Wtilde2.submat(sidx-1, sidx2-1, eidx-1, eidx2-1) = Wtilde;
    
  }
  
  arma::mat C2(J, n*J);
  C2.zeros();
  
  double xkr=0, xki=0;
  arma::rowvec ckr(J);
  ckr.zeros();
  
  
  for(int r=1;r<=J;r++){
    
    for(int k=1;k<=n;k++){
      xkr = X(k-1,r-1);
      ckr.zeros();
      
      tmp = X.row(k-1)*beta;
      
      
      
      for(int i=1;i<=J;i++){
        xki = X(k-1, i-1);
        ckr[i-1] = xki*xkr*fl_dot(tmp(0,0), strDistr=strDistr);
      }
      sidx= (k-1)*J+1;
      eidx = sidx+J-1;
      C2.submat(r-1, sidx-1, r-1, eidx-1) = ckr;
      
    }
  }
  

  
  arma::vec P2tilde(pow(n,2));
  P2tilde.zeros();
  arma::vec Pk(n);
  Pk.zeros();
  
  for(int k=1;k<=n;k++){
    
    sidx = (k-1)*n+1;
    eidx = sidx+n-1;
    Pk = P2.col(k-1);
    P2tilde.subvec(sidx-1, eidx-1) = Pk;
  }
  
  arma::vec d = C2*Wtilde2*P2tilde;
  return d;
  
  arma::vec out(J);
  return out;
  
}



double Get_AVec_r(arma::mat X, arma::mat D, arma::vec beta, int r, Rcpp::String strDistr="Logit"){
  
  int n = X.n_rows;
  int J = X.n_cols;
  
  arma::mat Wn = inv(Get_Wn(X, D, beta));
  
  arma::mat Dstar = D*D.t();
  
  arma::mat Pn(n,n);
  Pn.zeros();
  
  arma::mat Lambda_n(n,n);
  Lambda_n.zeros();
  
  arma::mat Lambda_dot_n(n,n);
  Lambda_dot_n.zeros();
  
  
  arma::mat tmp(1,1);
  
  
  arma::vec PnVec(n);
  arma::vec fnVec(n);
  
  PnVec.zeros();
  fnVec.zeros();
  
  double fk=0, Fk=0, fdotk=0, f2dotk=0;
  
  for(int k=1;k<=n;k++){
    
    tmp = X.row(k-1)*beta;
    fk = fl(tmp(0,0), strDistr=strDistr);
    Lambda_n(k-1,k-1) = fk;
    
    Lambda_dot_n(k-1,k-1) = fl_dot(tmp(0,0), strDistr=strDistr);
    
    Fk = Fl(tmp(0,0), strDistr=strDistr);
    
    Pn(k-1,k-1) = Fk*(1-Fk);
    PnVec[k-1] = Fk;
    fnVec[k-1] = fk;
    
  }
  
  arma::mat P2 = D*D.t()*Pn*D*D.t();
  
  arma::rowvec WVec(pow(J,2));
  
  int sidx=0, eidx=0;
  int sidx2=0, eidx2=0;
  
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*J+1;
    eidx = sidx+J-1;
    
    WVec.subvec(sidx-1, eidx-1) = Wn.row(i-1);
    
  }
  
  arma::mat LP2L = Lambda_n* P2* Lambda_n;
  
  arma::mat LP2L_Tilde(n*J, n*J);
  LP2L_Tilde.zeros();
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*n+1;
    eidx = sidx+n-1;
    
    for(int k=1;k<=J;k++){
      sidx2 = (k-1)*n+1;
      eidx2 = sidx2+n-1;
      
      LP2L_Tilde.submat(sidx-1, sidx2-1, eidx-1, eidx2-1) = LP2L;
    }
    
  }
  
  arma::mat Xtilde(n*J, J);
  Xtilde.zeros();
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*n+1;
    eidx = sidx+n-1;
    
    Xtilde.submat(sidx-1, i-1, eidx-1, i-1) = X.col(i-1);
  
  }
  
  arma::mat XLX = Xtilde.t()*LP2L_Tilde*Xtilde;
  
  arma::mat XLX2 = XLX;
  
  
  double c11 = 0;
  double c12 = 0;
  double c2 = 0;
  double c3 = 0;
  
  
  arma::mat xx_k(J,J);
  xx_k.zeros();
  
  arma::rowvec xk(J);
  xk.zeros();
  
  arma::mat XLX_k(pow(J,2), pow(J,2));
  
  arma::mat WXW(1,1);
  
  double tmpval=0, tmpval2=0;
  double xkr=0;
  
  
  
  arma::mat Pk_a(n,n);
  Pk_a.zeros();
  arma::mat Ak(n,n);
  Ak.zeros();
  
  
  arma::mat DLX(1,1);
  
  arma::mat Xbk(1,1);
  arma::mat Xbm(1,1);
  
  
  for(int k=1;k<=n;k++){
    Xbk = X.row(k-1)*beta;
    
    fdotk = fl_dot(Xbk(0,0), strDistr=strDistr);
    f2dotk = fl_2dot(Xbk(0,0), strDistr=strDistr);
    fk = fl(Xbk(0,0), strDistr=strDistr);
    
    
    
    DLX = Dstar.row(k-1)*Lambda_n*X.col(r-1);
    xk = X.row(k-1);
    xx_k = xk.t()* X.row(k-1);
    XLX_k = kron(XLX, xx_k);
    
    WXW = WVec*XLX_k*WVec.t();
    
    c11 -= DLX(0,0)*fdotk*WXW(0,0);
    
    
    ////////////// c12
    
    tmpval=0;
    tmpval2=0;
    
    xkr = X(k-1, r-1);
    
    
    for(int m=1;m<=n;m++){
      Xbm = X.row(m-1)*beta;
      
      double fdotm = fl_dot(Xbm(0,0), strDistr=strDistr);
      
      double fm = fl(Xbm(0,0), strDistr=strDistr);
      double Fm = Fl(Xbm(0,0), strDistr=strDistr);
      
      double dkm_star = Dstar(k-1,m-1);
      double xmr = X(m-1, r-1);
      
      Pk_a(m-1,m-1) = (Fm-3*pow(Fm,2)+2*pow(Fm, 3))*dkm_star;
      
      arma::mat xx_km = xk.t()*X.row(m-1);
      arma::mat XLX_km = kron(XLX, xx_km);
      WXW = WVec*XLX_km*WVec.t();
      tmpval += dkm_star*xmr*fdotm*WXW(0,0);
      
      
      /////// c2
      tmpval2 += dkm_star*fm*WXW(0,0);
      
    }
    
    c12 -= fk*tmpval;
    c2 -= fdotk*xkr*tmpval2;
    
    
    
    /////////////////c3
    
    Ak = Dstar*Pk_a*Dstar;
    
    arma::mat LAkL = Lambda_n*Ak*Lambda_n;
    arma::mat LAkL_Tilde(n*J, n*J);
    LAkL_Tilde.zeros();
    
    for(int i=1;i<=J;i++){
      
      sidx =(i-1)*n+1;
      eidx = sidx+n-1;
      
      for(int k2=1;k2<=J;k2++){
        sidx2 = (k2-1)*n+1;
        eidx2 = sidx2+n-1;
        
        LAkL_Tilde.submat(sidx-1, sidx2-1, eidx-1, eidx2-1)=LAkL;
      }
    }
    
    
    XLX2 = Xtilde.t()*LAkL_Tilde*Xtilde;
    
    XLX_k = kron(XLX2, xx_k);
    
    WXW = WVec*XLX_k*WVec.t();
    
    c3 += f2dotk*xkr*WXW(0,0);
    
  }
  

  double out = c11+c12+c2+c3;
  
  return out;
  
  
  
  
  
}


arma::vec Get_AVec(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr="Logit"){
  
  int n = X.n_rows;
  int J = X.n_cols;
  
  arma::mat Wn = inv(Get_Wn(X, D, beta));
  
  arma::mat Dstar = D*D.t();
  
  arma::mat Pn(n,n);
  Pn.zeros();
  
  arma::mat Lambda_n(n,n);
  Lambda_n.zeros();
  
  arma::mat Lambda_dot_n(n,n);
  Lambda_dot_n.zeros();
  
  
  arma::mat tmp(1,1);
  
  
  arma::vec PnVec(n);
  arma::vec fnVec(n);
  
  PnVec.zeros();
  fnVec.zeros();
  
  double fk=0, Fk=0, fdotk=0, f2dotk=0;
  
  for(int k=1;k<=n;k++){
    
    tmp = X.row(k-1)*beta;
    fk = fl(tmp(0,0), strDistr=strDistr);
    Lambda_n(k-1,k-1) = fk;
    
    Lambda_dot_n(k-1,k-1) = fl_dot(tmp(0,0), strDistr=strDistr);
    
    Fk = Fl(tmp(0,0), strDistr=strDistr);
    
    Pn(k-1,k-1) = Fk*(1-Fk);
    PnVec[k-1] = Fk;
    fnVec[k-1] = fk;
    
  }
  
  arma::mat P2 = D*D.t()*Pn*D*D.t();
  
  arma::rowvec WVec(pow(J,2));
  
  int sidx=0, eidx=0;
  int sidx2=0, eidx2=0;
  
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*J+1;
    eidx = sidx+J-1;
    
    WVec.subvec(sidx-1, eidx-1) = Wn.row(i-1);
    
  }
  
  arma::mat LP2L = Lambda_n* P2* Lambda_n;
  
  arma::mat LP2L_Tilde(n*J, n*J);
  LP2L_Tilde.zeros();
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*n+1;
    eidx = sidx+n-1;
    
    for(int k=1;k<=J;k++){
      sidx2 = (k-1)*n+1;
      eidx2 = sidx2+n-1;
      
      LP2L_Tilde.submat(sidx-1, sidx2-1, eidx-1, eidx2-1) = LP2L;
    }
    
  }
  
  arma::mat Xtilde(n*J, J);
  Xtilde.zeros();
  
  for(int i=1;i<=J;i++){
    
    sidx = (i-1)*n+1;
    eidx = sidx+n-1;
    
    Xtilde.submat(sidx-1, i-1, eidx-1, i-1) = X.col(i-1);
    
  }
  
  arma::mat XLX = Xtilde.t()*LP2L_Tilde*Xtilde;
  
  arma::mat XLX2 = XLX;
  
  
  double c11 = 0;
  double c12 = 0;
  double c2 = 0;
  double c3 = 0;
  
  
  arma::mat xx_k(J,J);
  xx_k.zeros();
  
  arma::rowvec xk(J);
  xk.zeros();
  
  arma::mat XLX_k(pow(J,2), pow(J,2));
  
  arma::mat WXW(1,1);
  
  double tmpval=0, tmpval2=0;
  double xkr=0;
  
  
  
  arma::mat Pk_a(n,n);
  Pk_a.zeros();
  arma::mat Ak(n,n);
  Ak.zeros();
  
  
  arma::mat DLX(1,1);
  
  arma::mat Xbk(1,1);
  arma::mat Xbm(1,1);
  
  
  arma::vec AVec(J);
  AVec.zeros();
  
  
  for(int r=1;r<=J;r++){
    
    
    c11 = 0;
    c12 = 0;
    c2 = 0;
    c3 = 0;
    
    
    for(int k=1;k<=n;k++){
      Xbk = X.row(k-1)*beta;
      
      fdotk = fl_dot(Xbk(0,0), strDistr=strDistr);
      f2dotk = fl_2dot(Xbk(0,0), strDistr=strDistr);
      fk = fl(Xbk(0,0), strDistr=strDistr);
      
      
      
      DLX = Dstar.row(k-1)*Lambda_n*X.col(r-1);
      xk = X.row(k-1);
      xx_k = xk.t()* X.row(k-1);
      XLX_k = kron(XLX, xx_k);
      
      WXW = WVec*XLX_k*WVec.t();
      
      
      c11 -= DLX(0,0)*fdotk*WXW(0,0);
      
      
      ////////////// c12
      
      tmpval=0;
      tmpval2=0;
      
      xkr = X(k-1, r-1);
      
      
      for(int m=1;m<=n;m++){
        Xbm = X.row(m-1)*beta;
        
        double fdotm = fl_dot(Xbm(0,0), strDistr=strDistr);
        
        double fm = fl(Xbm(0,0), strDistr=strDistr);
        double Fm = Fl(Xbm(0,0), strDistr=strDistr);
        
        double dkm_star = Dstar(k-1,m-1);
        double xmr = X(m-1, r-1);
        
        Pk_a(m-1,m-1) = (Fm-3*pow(Fm,2)+2*pow(Fm, 3))*dkm_star;
        
        arma::mat xx_km = xk.t()*X.row(m-1);
        arma::mat XLX_km = kron(XLX, xx_km);
        WXW = WVec*XLX_km*WVec.t();
        tmpval += dkm_star*xmr*fdotm*WXW(0,0);
        
        
        /////// c2
        tmpval2 += dkm_star*fm*WXW(0,0);
        
      }
      
      c12 -= fk*tmpval;
      c2 -= fdotk*xkr*tmpval2;
      
      
      
      /////////////////c3
      
      Ak = Dstar*Pk_a*Dstar;
      
      arma::mat LAkL = Lambda_n*Ak*Lambda_n;
      arma::mat LAkL_Tilde(n*J, n*J);
      LAkL_Tilde.zeros();
      
      for(int i=1;i<=J;i++){
        
        sidx =(i-1)*n+1;
        eidx = sidx+n-1;
        
        for(int k2=1;k2<=J;k2++){
          sidx2 = (k2-1)*n+1;
          eidx2 = sidx2+n-1;
          
          LAkL_Tilde.submat(sidx-1, sidx2-1, eidx-1, eidx2-1)=LAkL;
        }
      }
      
      
      XLX2 = Xtilde.t()*LAkL_Tilde*Xtilde;
      
      XLX_k = kron(XLX2, xx_k);
      
      WXW = WVec*XLX_k*WVec.t();
      
      c3 += f2dotk*xkr*WXW(0,0);
      
    }
    
    AVec[r-1] = c11+c12+c2+c3;
    
  }
  
  

  return AVec;
  
  
  
  
  
}


// [[Rcpp::export]]
arma::vec Get_bias(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr="Logit"){
  
  
  int n = X.n_rows;

  arma::mat Lambda_n(n,n);
  Lambda_n.zeros();
  
  double fk=0;
  
  arma::mat tmp(1,1);
  
  for(int k=1;k<=n;k++){
    tmp = X.row(k-1)*beta;
    
    fk = fl(tmp(0,0), strDistr=strDistr);
    Lambda_n(k-1,k-1) = fk;
  }
 
  arma::mat CMat = X.t()*Lambda_n*D*D.t()*Lambda_n*X;
  arma::vec dVec = Get_d(X, D, beta, strDistr=strDistr);
  
  
  arma::vec AVec = Get_AVec(X, D, beta, strDistr=strDistr);
  arma::vec biasVec = inv(CMat) * (dVec+0.5*AVec);
  

  return biasVec;
  
}




// [[Rcpp::export]]
List Find_MDBeta(arma::vec beta0, arma::vec Y, arma::mat X, arma::mat D, 
                     Rcpp::String strDistr, int nIter=100, double lr=0.01, double crit=1e-3,
                     bool bBias=false, bool bDisp=false){
  
  List lst = Get_beta(beta0, Y, X, D, strDistr, nIter=nIter, lr=lr, crit=crit, bDisp=bDisp);
  
  //if(bBias==true){
    
    //  arma::vec newbeta = lst[2];
    //arma::vec biasVec = Get_bias(X, D, newbeta, strDistr=strDistr);
    
    //newbeta -= biasVec;
    
    //lst[2] = newbeta;
  //}


  return lst;
  
}




