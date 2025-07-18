
#' Implementing a binary regression estimation.
#'
#' Performs L2-type optimization
#'@param formula  a description of the model to be fitted.
#'@param data a data frame that contains the response variable and predictors.
#'@param beta0 an initial value for beta to be estimated. The default is Null.. 
#'@param link a link function. It should be Logit, Probit, or LogLog.
#'@param bBias a logical value for the bias reduction. The default is FALSE
#'@param nIter the number of the iterations for the gradient decent method to find the MD estimator. The default value is 100.
#'@param lr the learning rate for the gradient descent method. The default value is 0.01.
#'@param crit the criterion used for exiting the iteration of the gradient descent method. The default value is 1e-3.
#'@param bDisp a logical value for displaying the iteration.  The default value is FALSE.

#'@return A list of the following values:
#'\describe{
#'\item{Iter_Num}{ the number of iterations taken to complete the estimation. If it is less than nIter, the convergence was made within the given iterations.}
#'
#'\item{diff}{the difference between the final estimators and the estimator of the previous stage, which shows the convergence status. If Iter_Num is less than nIter, then diff should be smaller than crit.}
#'
#'\item{beta_MDE}{the MD estimator.}
#'
#'\item{LossVec}{ a vector of values of the loss function at the estimator over the iterations.}
#'
#'\item{fitted_values}{ a vector of fitted values obtained using the link function and the MD estimator.}
#'
#'\item{residual}{ a vector of residuals, which is equal to the response substracted from fitted_values.}
#'
#'}
#'
#'@examples
#'####################

#'n=20
#'beta = c(1,-2,3.5)
#'p=3
#'
#'###### Generate a data frame including a response and design matrix
#'X = matrix(0, n,p)

#'for(i in 1:p){
#'  X[,i] = runif(n, -2, 2)
#'}

#'pVec = 1/( 1+exp(-X%*%beta) )
#'Y = rbinom(n, 1, pVec)

#'DM = cbind(Y,X)
#'colnames(DM) = c("Y", paste("X", 1:p, sep=""))
#'DM = data.frame(DM)
#'
#'link = "Logit"
#'lst = jwglm("Y~X1+X2+X3", data=DM, nIter=1000)
#'
#'@export
#'@importFrom Rcpp evalCpp
#'@useDynLib jwglm




jwglm = function(formula, data, beta0=FALSE, D=FALSE,  link="Logit", bBias=FALSE,
                 nIter=100, lr=0.01, crit=1e-3, bDisp=FALSE){
  
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  
  
  m <- match(c("formula", "data"),      names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  
  mt <- attr(mf, "terms") # allow model.frame to update it
  Y <- model.response(mf, "numeric")
  
  
  X <- model.matrix(mt, mf)
  
  XX = t(X)%*%X
  A = sqrtmat(XX, -0.5)
  
  
  if(is.null(ncol(D))==TRUE){
    
    D = X%*%A
    
  }
  
  if(length(beta0)==1){
    if(beta0==FALSE){
      Init_beta = rep(1, times=ncol(X))
      
    }else{
      Init_beta=beta0
    }
  }else{
    
    Init_beta = beta0
    
  }
  
  
  lst = Find_MDBeta(Init_beta, Y, X, D, strDistr=link, nIter=nIter, lr=lr, crit=crit,
                    bBias=bBias, bDisp=bDisp)
  
  
  dimm = dim(X)
  n=dimm[1]
  p = dimm[2]
  
  fittedVec = rep(0, times=n)
  
  beta_mde = lst[[3]]
  for(i in 1:n){
    xVec = matrix(X[i,], 1, p)
    xb = xVec%*%beta_mde
    fittedVec[i] = Fl(xb, strDistr=link)
  }
  
  residual = Y-fittedVec
  
  Final_list = list("Iter_Num" = lst[[1]], "diff"=lst[[2]], "beta_MDE" = lst[[3]], 
                    "LossVec"=lst[[4]], "fitted_values" = fittedVec, "residual" = residual)
  
  return(Final_list)
}