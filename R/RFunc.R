

sqrtmat<-function(MAT, EXP, tol=NULL){
  MAT <- as.matrix(MAT)
  matdim <- dim(MAT)
  if(is.null(tol)){
    tol=min(1e-7, .Machine$double.eps*max(matdim)*max(MAT))
  }
  if(matdim[1]>=matdim[2]){
    svd1 <- svd(MAT)
    keep <- which(svd1$d > tol)
    res <- t(svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep]))
  }
  if(matdim[1]<matdim[2]){
    svd1 <- svd(t(MAT))
    keep <- which(svd1$d > tol)
    res <- svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep])
  }
  return(res)
}

DGP = function(beta, n, p, bOut=FALSE){
  
  p = length(beta)
  
  X = matrix(0, n,p)
  X[,1] = rep(1, times=n)
  
  for(i in 2:p){
    X[,i] = runif(n, -2, 2)
    #X[,i] = rnorm(n, 1, 0.5)
  }
  
  pVec = 1/( 1+exp(-X%*%beta) )
  Y = rbinom(n, 1, pVec)
  
  if(bOut==TRUE){
    
    n2= floor(n*0.1)     ## original n2=5
    if(n2==0){
      n2=2
    }
    
    X2 = matrix(0, n2,p)
    X2[,1] = rep(1, times=n2)
    
    for(i in 2:p){
      
      X2[,i] = runif(n2, -2, 2)     ### original X2[,i] = runif(n, -7, -5)  
    }
    
    X2[,p] = runif(n2, 10, 20)
    
    #Y2 = rep(1, n2)
    
    pVec2 = 1/( 1+exp(-X2%*%beta) )
    Y2 = rbinom(n, 1, pVec2)
    
    
    X[(n-(n2-1)):n,] = X2
    Y[(n-(n2-1)):n] = Y2
  }
  
  #DM = data.frame(Y=Y, Intercept=X[,1], X1=X[,2], X2=X[,3])  
  #DM = data.frame("Y"=Y, "Intercept"=X[,1], "X1"=X[,2], "X2"=X[,3])
  
  tmpstr = "DM = data.frame(Y=Y, Intercept=X[,1]"
  
  for(i in 2:p){
    
    tmpstr=paste0(tmpstr, ", X", (i-1), "=X[,", i, "]" )
    
    if(i==p){
      tmpstr=paste0(tmpstr, ")" )
    }
  }
  eval(parse(text=tmpstr))
  
  lst = list(DM=DM, Y=Y, X=X)
  
  return(lst)  
}




