#' 
#' determine the number of change points by using LASSO
#' 
#' Author: Chuang Wan
#' 
#' Email: wanchuanghnu@126.com
#' 
#' @keywords SERseg; larsER
#' 
#' @param y A vector of response
#' @param x A scalar covariate with threshold
#' @param z A vector of covariates
#' @param tau the expectile level
#' @param tol  tolerance value, 1e-4 for default
#' @param max.iter the maximum iteration steps, 100 for default
#'
#' @return A list with the elements
#' \item{coef.est}{The estimated regression coefficients with intercept.}
#' \item{variance}{The estimated variance and covariance fot parmaters}
#' \item{iter}{The iteration steps.}
#'

#' @importFrom stats approxfun density lsfit sd
#' @importFrom Matrix nearPD
#' @importFrom lars
#'
#' @export
#'
#'
#'
require(lars)
library(lars)

#Iterative estimating procedures for change points

SERseg <- function(y,z,x,psi=NULL,k,grid,max.iter=100,tol=1E-4){
  
  n <- length(y)
  
  if(is.null(psi)) psi <- quantile(x,prob=seq(0,1,l=k+2)[-c(1,k+2)],names=FALSE)
  
  iter <- 1
  
  while(iter < max.iter){
    
    k <- length(psi)
    
    X <- matrix(rep(x,k),nrow=n)
    
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    
    U <- pmax((X-PSI),0)
    
    V <- ifelse((X>PSI),-1,0)

    XZ <- cbind(z,x,U,V)
    
    colnames(XZ) <- c("z","x",paste("U",1:ncol(U),sep=""),
                      paste("V",1:ncol(V),sep=""))
    
    obj<- expFit(XZ, y, tau,  max.iter=100, tol=1E-4)
    
    bet1 <- obj$coefficients
    
    #bet1.iter.se <- summary(obj)$coef[,2]
    bet1.iter.se <- diag(obj$variance)
    
    beta.U <- bet1[paste("U",1:ncol(U),sep="")]
    
    beta.U.se <- bet1.iter.se[paste("U",1:ncol(U),sep="")]
    
    beta.V <- bet1[paste("V",1:ncol(V),sep="")]
    
    beta.V.se <- bet1.iter.se[paste("V",1:ncol(V),sep="")]
    
    alpha <- bet1[c("Intercept","z","x")]
    
    alpha.se <- bet1.iter.se[c("Intercept","z","x")]
    
    psi.old <- psi
    
    psi <- psi.old + beta.V/beta.U
    
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    
    A <- apply((X <= PSI),2,all)
    
    B <- apply((X >= PSI),2,all)
    
    id.psi.ok <- !is.na((A+B)<=0)&(A+B)<=0
    
    #X <- X[,id.psi.ok,drop=FALSE]
    
    psi <- psi[id.psi.ok]
    
    bet.V <- beta.V[id.psi.ok]
    
    bet.V.se <- beta.V.se[id.psi.ok]
    
    bet.U <- beta.U[id.psi.ok]
    
    bet.U.se <- beta.U.se[id.psi.ok]
    ##################################  
    bet.V <- bet.V[order(psi)]
    
    bet.V.se <- bet.V.se[order(psi)]
    
    bet.U <- bet.U[order(psi)]
    
    bet.U.se <- bet.U.se[order(psi)]
    
    psi <- sort(psi)
    
    C <- (diff(psi)<grid)
    
    if(length(C[C==TRUE]) != 0){
      
      id.diff.no <-  which(C==TRUE)+1
      
      psi <- psi[-id.diff.no]
      
      bet.V <- bet.V[-id.diff.no]
      
      bet.V.se <- bet.V.se[-id.diff.no]
      
      bet.U <- bet.U[-id.diff.no]
      
      bet.U.se <- bet.U.se[-id.diff.no]
      
    }
    
    
    
    #########################
    
    #PSI <- PSI[,id.psi.ok,drop=FALSE]
    
    #U <- pmax((X-PSI),0)
    
    if(length(psi)<=0){
      
      obj$n.psi <- 0
      
      obj$psi <- NULL
      
      return(obj)
    }
    
    psi.se <- bet.V.se/abs(bet.U)
    
    iter <- iter +1
    
    if(max(abs(bet.V)) < tol) break
    
  }
  
  X <- matrix(rep(x,length(psi)),nrow=n)
  
  PSI <- matrix(rep(psi,rep(n,length(psi))),ncol=length(psi))
  
  U <- pmax((X-PSI),0)
  
  obj$psi <- sort(psi)
  
  obj$psi.se <- psi.se[order(psi)]
  
  beta.U <- bet.U[order(psi)]
  
  obj$est <- c(alpha,bet.U)
  
  obj$est.se <- c(alpha.se,bet.U.se)
  
  #obj$est.V <- bet.V[order(psi)]
  
  obj$U <- U
  
  obj$n.psi<- length(psi)
  
  obj$iter <- iter
  
  return(obj)
  
}


#Refine the candate change points by LASSO
LarsER <- function(y,z,x,k,psi,grid,max.iter=50,tol=1e-4){
  
  n <- length(y)
  
  if(length(x) != n) stop("Length of x and y differ")
  
  obj <- SERseg(y,z,x,psi,k,grid,max.iter=100,tol=1E-4)
  
  if(obj$n.psi == 0){
    
    ris <- list(n.psi=0,psi=NULL)
    
    #X <- cbind(1,x,z)
    
    #obj <- expFit(X, y, tau,  max.iter=100, tol=1E-4)
    
    return(ris)
    
  }
  
  psi0 <- obj$psi
  
  edf.psi <- TRUE
  
  type <- "bic"
  
  tipoAlg<- "lasso"
  
  Cn <- eval(parse(text="log(log(n))"))
  
  S <- 1
  
  olars <- lars(obj$U,y=y,type=tipoAlg,normalize=FALSE,intercept=TRUE,trace=FALSE)
  
  id.var.entry <- (1:ncol(obj$U))[order(olars$entry)]
  
  edf <- if(edf.psi) (olars$df-1)*2+1 else olars$df
  
  RSS <- olars$RSS
  
  min.r<-switch(type,
                bic = which.min(log(RSS/n)+log(n)*edf*Cn/n),
                #mdl = which.min(n*log(RSS/(n-edf))+Cn*pen.MDL(id.var.entry,as.numeric(table(-rowSums(obj$V))))),
                rss = max.rss(RSS)
  )
  crit<-switch(type,
               bic = (log(RSS/n)+log(n)*edf*Cn/n),
               #mdl = (n*log(RSS/(n-edf))+Cn*pen.MDL(id.var.entry,as.numeric(table(-rowSums(obj$V))))),
               rss = (RSS)
  )
  
  id <- sort(c(0,id.var.entry)[1:min.r])
  
  id <- id[-1]
  
  psi1 <- obj$psi[id]
  
  ris <- list(psi=psi1,n.psi=length(psi1))
  
  return(ris)
}


