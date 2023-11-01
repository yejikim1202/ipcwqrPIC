#' @importFrom stats as.formula binomial predict sd
NULL
#' Fit the doubly-censored AFT model with quantile regression model
#'
#' Fit inverse weighted quantile regression with doubly interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param T exactly observed time.
#' @param delta censoring indicator, 1: exactly observed; 2: right-censored; 3: left-censored;
#' @param x X matrix of baseline covariates.
#' @param tau quantile level.
#' @param estimation estimating method of partly interval censored, if estimation="DR", doubly robust estimator is estimated.
#' @param var.estimation variance estimating method, if var.estimation="IS", the induced smoothing method is used, and else if var.estimation="Bootstrap", variance bootstrapping method is used.
#' @param wttype weight estimating method, default is "KM", Beran's nonparametric KM estimating method as "Beran", and  Ishwaran's random survival forests KM estimating method as "Ishwaran".
#' @param hlimit bandwidth value, default is NULL.
#' @param contx1.pos position of the continuous covariate of variable x used in the kernel of the Beran method, the default is 1.
#' @param contx2.pos position of the continuous covariate of variable x used in the kernel of the Beran method, the default is 1.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param index index of cluster weight, default is 1
#' @param B the number of iterations in the bootstrap method., default is 100
#' @param maxit maximum time value of event time T or L and R, default is 100.
#' @param max.iter maximum number of iteration for the quantile regression estimator, default is 100.
#' @param tol tolerance of iteration for the quantile regression estimator, default is 1e-3.
#'
#' @return \code{dcrq} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{tau}: quantile level.
#'   \item \code{coefficients}: regression estimator.
#'   \item \code{se}: standard error estimates for \code{est}.
#'   \item \code{pvalue}: p-value.
#'   \item \code{lower bd}: lower bound of coefficients under 95% confidence level.
#'   \item \code{upper bd}: upper bound of coefficients under 95% confidence level.
#' }
#'
#' @details
#' see Kim et al., (2023+) for detailed method explanation.
#'
#' @references
#' 
#' Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.
#' 
#' Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.
#' 
#' Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time- dependent covariates. Biostatistics. 18 (1), 132–146.
#' 
#' Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.
#' 
#' Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#' 
#' Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). Inverse weighted quantile regression with partially interval-censored data.
#' 
#'
#' @examples
#' \dontrun{
#' # Simulations
#' n=200; x1=runif(n,-1.2,1.7); x2=rbinom(n,1,0.6)
#' T = 1.7+x1+x2+rnorm(n)*(1-0.1*x2)
#' L=runif(n,-2.8,1.9); R=L+runif(n,4.2,8.1)
#' Y=pmin(R,pmax(T,L))
#' delta=case_when(
#'  T<L ~ 3,
#'  T>R ~ 2,
#'  TRUE ~ 1 #observed
#')
#'L=L; R=R; T=T; delta=delta; x=cbind(x1,x2); tau=0.3
#'dcrq(L,R,T,delta,x,tau,var.estimation="IS")
#'dcrq(L,R,T,delta,x,tau,var.estimation="Bootstrap")
#'dcrq(L,R,T,delta,x,tau,estimation = "DR",var.estimation="IS")
#'dcrq(L,R,T,delta,x,tau,wttype = "Ishwaran",var.estimation="IS")
#'dcrq(L,R,T,delta,x,tau,wttype = "Beran",hlimit = 0.1,var.estimation="IS")
#' }
#' @export
#'
#'
#'


dcrq=function(L,R,T,delta,x,tau,estimation=NULL,var.estimation=NULL,wttype="KM",hlimit=NULL,contx1.pos=1,contx2.pos=1,id=NULL,index=1,B=100,maxit=100,max.iter=100,tol=1e-3){
  
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(quantreg)
  library(glmnet)
  library(randomForestSRC)
  
  
  wtft = function(L,R,T=NULL,estimation=NULL,delta){
    
    if(sum(delta==4)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
      deltaL = ifelse(delta==4|delta==3,1,0)
      deltaR = ifelse(delta==4|delta==2,1,0)
      kml = survfit(Surv(-L,deltaL==0)~1)
      kmr = survfit(Surv(R,deltaR==0)~1)
      
    }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
      deltaL = ifelse(delta==3,1,0)
      deltaR = ifelse(delta==2,1,0)
      kml = survfit(Surv(-Y,deltaL==1)~1)
      kmr = survfit(Surv(Y,deltaR==1)~1)
      
    }else if(sum(delta==2)!=0 & sum(delta==3)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y)
      deltaR = ifelse(delta==2,1,0)
      kmr = survfit(Surv(Y,deltaR==1)~1)
      
    }else if(sum(delta==3)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y)
      deltaL = ifelse(delta==3,1,0)
      kml = survfit(Surv(-Y,deltaL==1)~1)
    }
    
    ww = rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        
        if(sum(delta==1)==n){
          ww[i] = 1
        }else if(sum(delta==4)!=0){ 
          sl=approx(c(0,-kml$time,100),c(1,1-kml$surv,0), xout=L[i])$y
          sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
          ww[i] = 1/pmax(1-(sr-sl), 1e-3)
        }else if(sum(delta==2)!=0 & sum(delta==3)!=0 & is.null(estimation)){
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr-sl, 1e-3)
        }else if(sum(delta==2)!=0 & sum(delta==3)!=0 & estimation=="DR"){ 
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr, 1e-3)
        }else if(sum(delta==2)!=0 & sum(delta==3)==0){  
          sr = approx( c(0, kmr$time, maxit), c(1, kmr$surv, 0), xout=Y[i])$y
          ww[i] = 1/pmax(sr, 1e-3)
        }else if(sum(delta==3)!=0){ 
          sl = approx( c(0, -kml$time, maxit), c(1, 1-kml$surv,0), xout=Y[i])$y
          ww[i] = 1/pmax(1-sl, 1e-3)
        }
      }
    }
    ww
  }
  
  Berfunc = function(L,R,T=NULL,x,delta) {
    
    if(sum(delta==1)==n){
      Y=T; n=length(Y); y=Y
      
    }else if(sum(delta==4)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); Y=pmax(ifelse(delta==4,R,L),1e-8); n=length(L); y=Y
      
    }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y); y=Y
      deltaL = ifelse(delta==3,1,0)
      deltaR = ifelse(delta==2,1,0)
      
    }else if(sum(delta==2)!=0 & sum(delta==3)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y); y=Y
      deltaR = ifelse(delta==2,1,0)
      
    }else if(sum(delta==3)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y); y=Y
      deltaL = ifelse(delta==3,1,0)
      
    }
    
    ker = dnorm(outer(x[,contx1.pos],x[,contx2.pos],"-")/hlimit)
    Wnj = ker / rowSums(ker)
    if(sum(delta==4)!=0){
      denomr = rowSums(outer(y,y,">=")*(Wnj))
      denoml = rowSums(outer(y,y,"<=")*(Wnj))
    }else{
      denomr = rowSums(outer(y,y,"<=")*(Wnj))
      denoml = rowSums(outer(y,y,">=")*(Wnj))
    }
    
    sr = sl = ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        
        y0 = y[i]; nom = Wnj[,i]
        
        if(sum(delta==1)==n){
          ww[i] = 1
        }else if(sum(delta==4)!=0){ 
          etar = 1*(y<=y0 & delta!=1)
          etal = 1*(y>=y0 & delta!=1)
          sr = prod((1 - nom/denomr)^etar)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax(1-(sr-sl), 1e-3)
          
        }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
          etar = 1*(y<=y0 & deltaR==1)
          etal = 1*(y>=y0 & deltaL==1)
          sr = prod((1 - nom/denomr)^etar)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax( sr-sl, 1e-3)
          
        }else if(sum(delta==2)!=0 & sum(delta==3)==0){  
          etar = 1*(y<=y0 & deltaR==1)
          sr = prod((1 - nom/denomr)^etar)
          ww[i] = 1/pmax(sr, 1e-3)
          
        }else if(sum(delta==3)!=0){ 
          etal = 1*(y>=y0 & deltaL==1)
          sl = 1-prod((1 - nom/denoml)^etal)
          ww[i] = 1/pmax(1-sl, 1e-3)
          
        }
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  
  
  Ishfunc = function(L,R,T=NULL,x,delta) {
    
    if(sum(delta==4)!=0){ 
      #pic
      L = pmax(L,1e-8); R = pmax(R,1e-8); n=length(L)
      deltaL = ifelse(delta==4|delta==3,0,1)
      deltaR = ifelse(delta==4|delta==2,0,1)
      dt=data.frame(L=L,R=R,statusl=deltaL,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
      #dc
      L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
      deltaL = ifelse(delta==3,0,1)
      deltaR = ifelse(delta==2,0,1)
      dt=data.frame(L=L,R=R,statusl=deltaL,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==2)!=0 & sum(delta==3)==0){  
      #rc
      R=pmax(R,1e-8); Y=pmin(R,T); n=length(Y); y=Y
      deltaR = ifelse(delta==2,0,1)
      dt=data.frame(R=R,statusr=deltaR,x=x,xx=1)
      
    }else if(sum(delta==3)!=0){ 
      #lc
      L = pmax(L,1e-8); Y=pmax(L,T); n=length(Y); y=Y
      deltaL = ifelse(delta==3,0,1)
      dt=data.frame(L=L,statusl=deltaL,x=x,xx=1)
    }
    
    
    if(sum(delta==0)==n){
      survl=survr=0
      
    }else if(sum(delta==4)!=0){ 
      kml.obj <- rfsrc(Surv(-L, statusl==1) ~ .-L-statusl-R-statusr, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
      
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-L-statusl-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
      kml.obj <- rfsrc(Surv(-L, statusl==1) ~ .-L-statusl-R-statusr, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
      
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-L-statusl-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==2)!=0 & sum(delta==3)==0){
      kmr.obj <- rfsrc(Surv(R, statusr==0) ~ .-R-statusr, data=dt)
      # kmr.obj <- predict(rfsrc(Surv(R, statusr) ~ xx, data=dt))
      kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
      survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
      
    }else if(sum(delta==3)!=0){ 
      kml.obj <- rfsrc(Surv(L, statusl==0) ~ .-L-statusl, data=dt)
      # kml.obj <- predict(rfsrc(Surv(L, statusl) ~ xx, data=dt))
      kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
      survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    }
    
    
    ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        
        if(sum(delta==1)==n){
          ww[i] = 1
          
        }else if(sum(delta==4)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=L[i])$y
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=R[i])$y
          ww[i] = 1/pmax(1-(sr-sl),1e-3)
          
        }else if(sum(delta==2)!=0 & sum(delta==3)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=Y[i])$y
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=Y[i])$y
          ww[i] = 1/pmax( sr-sl, 1e-3)
          
        }else if(sum(delta==2)!=0 & sum(delta==3)==0){ 
          sr = approx( c(0, (kmr$event.info$time.interest), maxit), c(1, survr, 0), xout=Y[i])$y
          ww[i] = 1/pmax(sr, 1e-3)
          
        }else if(sum(delta==3)!=0){ 
          sl = approx( c(0, (kml$event.info$time.interest), maxit), c(1, survl, 0), xout=Y[i])$y
          ww[i] = 1/pmax(1-sl, 1e-3)
          
        }
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  DCrq=function(L,R,T,x,delta,ww,tau,eta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
    rq((Y)~x, weights = ww*eta, tau = tau)$coef #int, beta1, beta2
  }
  
  
  Efunc=function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) )
    U = as.vector( t(xx *eta)%*%(Phi* ww  - tau) )
    U/cluster
  }
  
  DREfunc=function(L,R,T,x,delta,tau,ww,wr,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
    wl=(ww-wr)/(ww*wr); wl[is.nan(wl)]=0; n=length(Y); 
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss =  pmax(1e-3, sqrt(diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wwind = ww*ind
    U = as.vector( t(xx *eta)%*%(wwind - tau) )
    UR=matrix(0,p,1)
    UL=matrix(0,p,1)
    
    for (i in 1:n) {
      yind=Y>=Y[i]
      denom=sum(yind*eta ) 
      if(delta[i]==2){
        indr=Y>=Y[i]
        dNir = Y<=Y[i]
        resr = as.numeric((Y - xx%*%beta)*indr)
        ind2 = ifelse(resr<=0,1,0)
        dMr=dNir-( (yind/denom) *dNir)
        Rft=(t(xx*wr*dMr*indr*eta)%*%( ind2 - tau ))
        UR=UR+((Rft/n))
      }
      if(delta[i]==3){
        indl=Y<=Y[i]
        dNil = Y>=Y[i]
        resl = as.numeric((Y - xx%*%beta)*indl)
        ind3 = ifelse(resl<=0,1,0)
        dMl=dNil-( ((1-yind)/(n+1-denom)) *dNil)
        Lft=(t(xx*wl*dMl*indl *eta)%*%( ind3 - tau ))
        UL=UL+((Lft/n))
      }
    }
    (U+(UR)+(UL))/(cluster)
  }
  
  
  Afunc=function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    exactnum=sum(ifelse(delta==1,1,0))
    ss =  sqrt(pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    wwss=ww/ss
    phi = as.vector( dnorm( -res/ss )* (wwss))
    A = t(phi * xx  *eta) %*% xx + diag(p)*0.05
    A/cluster
  }
  
  Gfunc=function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) )
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wwind = ww*ind
    
    Gam=(xx)*as.numeric(wwind-tau)
    Gamma=( t(Gam)%*%(Gam*eta) )
    GammaR=matrix(0,p,p)
    GammaL=matrix(0,p,p)
    
    for(i in 1:n){
      yind=(Y>=Y[i])
      denom=sum(yind*eta)
      if(delta[i]==2){
        dNir = (Y<=Y[i])*eta
        Br=t(xx*eta)%*%(yind*wwind)
        Br2=t(xx*eta)%*%(yind*wwind*dNir)
        R = (Br)/(denom)
        R2 = (Br2)/(denom)
        GammaR=GammaR+(R)%*%t(R2)
      }
      if(delta[i]==3){
        dNil = (Y>=Y[i])*eta
        Bl=t(xx*eta)%*%((1-yind)*wwind)
        Bl2=t(xx*eta)%*%((1-yind)*wwind*dNil)
        L = (Bl)/((n+1-denom))
        L2 = (Bl2)/((n+1-denom))
        GammaL=GammaL+(L)%*%t(L2)
      }
    }
    (Gamma-GammaR+GammaL)/cluster
  }
  
  Gfunc2 = function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma) {
    n=length(L);
    library(MASS)
    Shat = t(replicate(B,{
      id = sample(n,n,replace = TRUE)
      Efunc(L=L[id],R=R[id],T=T[id],x=x[id,],delta=delta,tau=tau,ww=ww[id],eta=eta[id],cluster=cluster,beta = beta, Sigma = Sigma)
    }))
    Var = cov(Shat) * (cluster)
    Var
  }
  
  Gfunc3= function(L,R,T,x,delta,tau,ww,eta,id,cluster,beta,Sigma) {
    n=length(L)
    library(MASS)
    Shat = t(replicate(B,{
      tabid=as.vector(table(id))
      idx = as.vector(unlist(lapply(tabid, function(x) sample(x=x,size=x,replace = TRUE))))
      Efunc(L=L[idx],R=R[idx],T=T[idx],x=x[idx,],delta=delta,tau=tau,ww=ww[idx],eta=eta[idx],cluster=cluster,beta = beta, Sigma = Sigma)
    }))
    Var = cov(Shat) * (cluster)
    Var
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc,cluster){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma/n
  }
  
  L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmin(R,T) );n=length(Y);
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else if(index==0){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  
  if(wttype=="KM"){ww=wtft(L=L,R=R,T=T,delta=delta);}
  if(wttype=="Ishwaran"){ww=Ishfunc(L=L,R=R,T=T,delta=delta,x=x);}
  if(wttype=="Beran"){ww=Berfunc(L=L,R=R,T=T,delta=delta,x=x);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = DCrq(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,tau=tau)
  old_Sigma = Sigma = diag(p)/cluster
  
  
  i=0; eps=1;
  while (i<max.iter & eps >= tol ) {
    Amat = Afunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    if(is.null(estimation)){
      new_beta = c(old_beta) - solve(Amat)%*%Efunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/n
    }else if(estimation=="DR"){
      wr=wtft(L=L,R=R,T=T,estimation="DR",delta=delta)# wr=Rwtfunc(L=L,R=R,T=T,delta=delta)
      new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,wr=wr,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/n
    }
    
    if(var.estimation=="IS"){
      Gamma = Gfunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }else if(var.estimation=="Bootstrap" & is.null(id)){
      Gamma = Gfunc2(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }else if(var.estimation=="Bootstrap" & is.null(id)==F){
      Gamma = Gfunc3(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,id=id,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }
    new_Sigma = up_Sigma(Y=Y,Afunc=Amat,Gfunc=Gamma,cluster=cluster)
    
    if (det(new_Sigma) <= 0) {
      new_beta = old_beta; new_Sigma = old_Sigma
    }
    eps = max(max(abs(new_beta - old_beta)),
              max(abs(new_Sigma - old_Sigma)))
    old_beta = new_beta; old_Sigma = new_Sigma; 
    i = i+1
  }
  se = sqrt(diag(new_Sigma))
  res=data.frame(tau=c(rep(tau,p)),
                 est=new_beta,se=se,
                 pvalue = 1 - pnorm(abs(new_beta/se)),
                 lb = new_beta-1.96*se, ub = new_beta+1.96*se)
  colnames(res)=c("tau","coefficients","se","pvalue","lower bd","upper bd")
  rownames(res)[1]="Intercept"
  round((res), 6) 
}

Etable=function(est,se,beta0){
  bias=colMeans(est)-beta0
  ase=apply(est, 2, sd)
  ese=colMeans(se)
  cpl=t( t(est)-1.96*t(se) < beta0 )
  cpu=t( t(est)+1.96*t(se) > beta0 )
  cp=colMeans( cpl*cpu )
  (round(cbind(bias,ase,ese,cp), 3))
}
