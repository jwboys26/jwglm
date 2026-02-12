#ifndef BIAS_H
#define BIAS_H


arma::vec Get_d(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr);
double Get_AVec_r(arma::mat X, arma::mat D, arma::vec beta, int r, Rcpp::String strDistr);
arma::vec Get_AVec(arma::mat X, arma::mat D, arma::vec beta, Rcpp::String strDistr);



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


#endif
