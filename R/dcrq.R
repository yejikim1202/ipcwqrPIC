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
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param index index of cluster weight, default is 1
#' @param maxit maximum number of iteration for the log-rank estimator, default is 100.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-3.
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

dcrq=function(L,R,T,delta,x,tau,estimation=NULL,var.estimation=NULL,wttype="KM",hlimit=NULL,id=NULL,index=1,maxit=100,tol=1e-3){
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(extRemes)
  library(quantreg)
  library(randomForestSRC)
  
  wtfunc=function(L,R,T,delta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    ww
  }
  
  
  Rwtfunc=function(L,R,T,delta){
    
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
    n=length(Y);
    kml = survfit(Surv(-Y,delta==3)~1)
    kmr = survfit(Surv(Y,delta==2)~1)
    wr = rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, -kml$time, 100), c(1,1-kml$surv,0), xout=Y[i])$y
        sr = approx( c(0, kmr$time, 100), c(1, kmr$surv, 0), xout=Y[i])$y
        wr[i] = 1/pmax( sr, 1e-3)
      }
    }
    wr
  }
  
  
  Berwtfunc = function(L,R,T,x,delta, h=NULL) {
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y); y=Y
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,"<=")*(Wnj))
    denoml = rowSums(outer(y,y,">=")*(Wnj))
    for (i in 1:n) {
      if(delta[i]==1){
        y0 = y[i]
        etar = 1*(y<=y0 & delta==2)
        etal = 1*(y>=y0 & delta==3)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    srl
  }
  
  Ishrfwtfunc = function(L,R,T,x,delta) {
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    status=ifelse(delta==1,1,0)
    dt=data.frame(L=L,R=R,status=status)
    
    kml.obj <- rfsrc(Surv(L, status) ~ ., data=dt)
    kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
    survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    
    kmr.obj <- rfsrc(Surv(R, status) ~ ., data=dt)
    kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
    survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
    
    ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, (kml$time), 100), c(1,survl,0), xout=Y[i])$y
        sr = approx( c(0, (kmr$time), 100), c(1, survr, 0), xout=Y[i])$y
        ww[i] = 1/pmax( sr-sl, 1e-3)
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  DCrq=function(L,R,T,x,delta,ww,tau,eta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    rq((Y)~x, weights = ww*eta, tau = tau)$coef #int, beta1, beta2
  }
  
  
  Efunc=function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss = sqrt( pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) )
    U = as.vector( t(xx *eta)%*%(Phi* ww  - tau) )
    U/cluster
  }
  
  Efunc2=function(L,R,T,x,delta,tau,ww,eta,cluster,beta){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wwind = ww*ind
    U = as.vector( t(xx *eta*ww)%*%(ind  - tau) )
    U/cluster
  }
  
  DREfunc=function(L,R,T,x,delta,tau,ww,wr,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    wl=(ww-wr)/(ww*wr); wl[is.nan(wl)]=0; wl[is.infinite(wl)]=0; n=length(Y); 
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
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y)
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss =  sqrt(pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    wwss=ww/ss
    phi = as.vector( dnorm( -res/ss )* (wwss))
    A = t(phi * xx  *eta) %*% xx + diag(p)*0.05
    A/cluster
  }
  
  Gfunc=function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
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
  
  Gfunc2 = function(L,R,T,x,delta,tau,ww,eta,cluster,beta,Sigma,B=100) {
    n=length(L);
    library(MASS)
    Shat = t(replicate(B,{
      id = sample(n,n,replace = TRUE)
      Efunc(L=L[id],R=R[id],T=T[id],x=x[id,],delta=delta,tau=tau,ww=ww[id],eta=eta[id],cluster=cluster,beta = beta, Sigma = Sigma)
    }))
    Var = cov(Shat*n)
    solve(Var)
  }
  
  
  Gfunc3= function(L,R,T,x,delta,tau,ww,eta,id,cluster,beta,Sigma,B=100) {
    n=length(L)
    # tabid=as.vector(table(id))
    # idx = as.vector(unlist(lapply(tabid, function(x) sample(x=x,size=x,replace = TRUE))))
    df=data.frame(num=1:n, id=id, L=L, R=R, T=T, ww=ww, x=x, eta=eta, delta=delta)
    len=length(unique(df$id))
    uniqid=unique(df$id)
    Shat2=matrix(0,ncol = 3,nrow = 1)
    for (i in 1:len) {
      df2 = subset(df,id==uniqid[i])
      library(MASS)
      Shat = t(replicate(B,{
        ndf2=nrow(df2)
        newx=as.matrix(cbind(df2$x.x1,df2$x.x2))
        idx = sample(ndf2,ndf2,replace = TRUE)
        Efunc(L=df2$L[idx],R=df2$R[idx],T=df2$T[idx],x=newx[idx,],delta=df2$delta[idx],tau=tau,ww=df2$ww[idx],eta=df2$eta[idx],cluster=cluster,beta = beta, Sigma = Sigma)
      }))
      Shat2=rbind(Shat2,(Shat))
    }
    Var = cov(Shat2*n)
    solve(Var)/sqrt(cluster)
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc,cluster){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma/n
  }
  
  L = pmax(L,1e-8); R=pmax(R,1e-8); Y=ifelse(L<R, pmin(R,pmax(L,T)), pmax(R,pmax(L,T)) );n=length(Y);
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  
  if(wttype=="Param"){ww=wtfunc(L=L,R=R,T=T,delta=delta);}
  else if(wttype=="Ishwaran"){ww=Ishrfwtfunc(L=L,R=R,T=T,delta=delta,x=x);}
  else if(wttype=="Beran" & is.null(hlimit)==F){ww=Berwtfunc(L=L,R=R,T=T,delta=delta,x=x,h=hlimit);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = DCrq(L=L,R=R,T=T,delta=delta,x=x,ww=ww,eta=eta,tau=tau)
  old_Sigma = Sigma = diag(p)/cluster
  
  
  i=0; eps=1; max.iter=100; tol = 1e-3; 
  while (i<max.iter & eps >= tol ) {
    Amat = Afunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    if(is.null(estimation)){
      new_beta = c(old_beta) - solve(Amat)%*%Efunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/n
    }
    else if(estimation=="DR"){
      wr=Rwtfunc(L=L,R=R,T=T,delta=delta)
      new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,wr=wr,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/n
    }
    
    if(var.estimation=="IS"){
      Gamma = Gfunc(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
      new_Sigma = up_Sigma(Y=Y,Afunc=Amat,Gfunc=Gamma,cluster=cluster)
    }
    else if(var.estimation=="Bootstrap" & is.null(id)){
      # new_beta = BB::dfsane(par=beta,fn=Efunc2,L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,control=list(trace=FALSE))$par
      new_Sigma = Gfunc2(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }
    else if(var.estimation=="Bootstrap" & is.null(id)==F){
      # new_beta = BB::dfsane(par=beta,fn=Efunc2,L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,control=list(trace=FALSE))$par
      new_Sigma = Gfunc3(L=L,R=R,T=T,x=x,delta=delta,tau=tau,ww=ww,eta=eta,id=id,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }
    
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
  colnames(res)=c("tau","coefficients","se","pvalue","95% lower bd","95% upper bd")
  rownames(res)[1]="Intercept"
  round((res), 6) 
}
