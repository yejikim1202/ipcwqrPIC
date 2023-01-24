#' @importFrom stats as.formula binomial lm predict sd
NULL
#' Fit the partly interval-censored AFT model with quantile regressions
#'
#' Fit inverse weighted quantile regression with partially interval-censored data
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param delta censoring indicator, 1: observed; 0: interval-censored.
#' @param x X matrix of baseline covariates.
#' @param tau quantile level.
#' @param estimation estimating method of partly interval censored, default is "NULL" and if estimation="dr", doubly robust estimator is estimated.
#' @param var.estimation estimating method of variance, if var.estimation="IS", induced smoothing method is used, else if var.estimation="bootstrap", bootstrap estimating method of Zeng and Lin is used.
#' @param wttype weight estimating method, default is "param", Beran's nonparametric KM estimating method as "Beran", and Ishwaran's nonparametric survival random forests estimating method as "Ishwaran"
#' @param hlimit bandwidth value, set \code{hlimit=NULL}.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param index index of cluster weight.
#' @param maxit maximum number of iteration for the IPCW method estimator, default is 100.
#' @param tol tolerance of iteration for the IPCW method estimator, default is 1e-3.
#'
#' @return \code{picrq} returns a data frame containing at least the following components:
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
#' see Kim et al., (2022+) for detailed method explanation.
#'
#' @references
#' Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.
#' 
#' Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.
#' 
#' Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time-dependent covariates. Biostatistics. 18 (1), 132–146.
#' 
#' Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.
#' 
#' Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.
#' 
#' Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). Inverse weighted quantile regression with partially interval-censored data.
#' 
#' Pan, C. (2021). PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. https://CRAN.R-project.org/package=PICBayes.
#'
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = runif(n,-1,1)
#' x2 = rbinom(n,1,0.43)
#' x = cbind(x1,x2)
#' T = 2 + x1 + x2 + rnorm(n)
#' U = (1 - 0.25*x1)*runif(n, -6, 5)
#' V = U + (1 - 0.1*x2)*runif(n, 6, 20)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' delta = ifelse(U==V, 1, 0)
#' tau=0.3
#' picrq(L=V,R=U,delta=delta,x=x,tau=tau,wttype="param",var.estimation="IS")
#' picrq(L=V,R=U,delta=delta,x=x,tau=tau,wttype = "Beran",hlimit=0.1,var.estimation="IS")
#' picrq(L=V,R=U,delta=delta,x=x,tau=tau,wttype = "Ishwaran",var.estimation="IS")
#' 
#' # Data example
#' library(PICBayes)
#' data("mCRC")
#' d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
#'                                       V = ifelse(y==2,L,R),
#'                                       # Cluster weighted data
#'                                       id=(rep(c(table(SITE)),c(table(SITE)))),
#'                                       # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
#'                                       x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
#'                                                     TRT_C == 1 ~ 1),
#'                                       # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
#'                                       x2= case_when(KRAS_C == 0 ~ 1,
#'                                                     KRAS_C == 1 ~ 0),
#'                                       delta = case_when(IC == 0 ~ 1,
#'                                                         IC == 1 ~ 0)
#'));
#'L=d$U;R=d$V; delta=d$delta
#'L=(log(d$U));R=log(d$V); delta=d$delta
#'x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
#'picrq(L,R,delta,x=x,tau=tau,id=id,index=1,wttype="param",var.estimation="IS")
#'picrq(L,R,delta,x=x,tau=tau,id=id,index=1,wttype="param",var.estimation="bootstrap")
#'picrq(L,R,delta,x=x,tau=tau,id=id,index=1,wttype = "Beran",hlimit = 0.9,var.estimation="bootstrap")
#' }
#' @export
#'
#'
library(extRemes)

picrq=function(L,R,delta,x,tau,estimation=NULL,var.estimation=NULL,wttype="param",hlimit=NULL,id=NULL,index=1,maxit=100,tol=1e-3){
  
  library(extRemes)
  library(MASS)
  library(tidyverse)
  library(survival)
  library(tidyverse)
  library(extRemes)
  library(quantreg)
  
  
  wtpicft=function(L,R,delta){
    
    L = pmax(L,1e-8); n=length(L)
    kml = survfit(Surv(L) ~ 1)
    kmr = survfit(Surv(R) ~ 1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
        ww[i] = 1/pmax((1-(sr-sl)),0.001)
      }
    }
    ww
  }
  
  
  Rwtpicft=function(L,R,delta){
    
    L = pmax(L,1e-8); n=length(L)
    kml = survfit(Surv(L) ~ 1)
    kmr = survfit(Surv(R) ~ 1)
    ww = rep(0,n)
    
    for (i in 1:n) {
      if(delta[i]==1){
        sl=approx(c(0,kml$time,100),c(1,kml$surv,0), xout=L[i])$y
        sr=approx(c(0,kmr$time,100),c(1,kmr$surv,0), xout=R[i])$y
        ww[i] = 1/pmax((sr),0.001)
      }
    }
    ww
  }
  
  Berwtpicfunc = function(L,R,x,delta, h=NULL) {
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); y=Y
    ker = dnorm(outer(x[,1],x[,1],"-")/h)
    Wnj = ker / rowSums(ker)
    sr = sl = srl= rep(0,n)
    denomr = rowSums(outer(y,y,"<=")*(Wnj))
    denoml = rowSums(outer(y,y,">=")*(Wnj))
    for (i in 1:n) {
      if(delta[i]==1){
        y0 = y[i]
        etar = 1*(y<=y0 & delta==0)
        etal = 1*(y>=y0 & delta==0)
        nom = Wnj[,i]
        sr = prod((1 - nom/denomr)^etar)
        sl = 1-prod((1 - nom/denoml)^etal)
        srl[i] = 1/pmax(1-(sr-sl),0.001)
      }
    }
    srl[is.na(srl)]=0
    srl
  }
  
  Ishrfwtpicfunc = function(L,R,x,delta) {
    library(randomForestSRC)
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y);
    statusl=ifelse(delta==0,0,1)
    statusr=ifelse(delta==0,0,1)
    dt=data.frame(L=L,R=R,statusl=statusl,statusr=statusr)
    
    kml.obj <- rfsrc(Surv(L, statusl) ~ ., data=dt)
    kml <- get.brier.survival(kml.obj, cens.model="rfsrc")
    survl=kml$surv.aalen; survl[is.na(survl)]=0; survl
    
    kmr.obj <- rfsrc(Surv(R, statusr) ~ ., data=dt)
    kmr <- get.brier.survival(kmr.obj, cens.model="rfsrc")
    survr=kmr$surv.aalen; survr[is.na(survr)]=0; survr
    
    ww= rep(0,n)
    for (i in 1:n) {
      if(delta[i]==1){
        sl = approx( c(0, (kml$time), 100), c(1,survl,0), xout=Y[i])$y
        sr = approx( c(0, (kmr$time), 100), c(1, survr, 0), xout=Y[i])$y
        ww[i] = 1/pmax((1-(sr-sl)),0.001)
      }
    }
    ww[is.na(ww)]=0
    ww
  }
  
  PICrq=function(L,R,x,delta,tau,ww,eta){
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    rq((Y)~x, weights = ww*eta, tau = tau)$coef #int, beta1, beta2
  }
  
  Efunc=function(L,R,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y);
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss =  sqrt(pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    Phi = as.vector( pnorm( -res/ss ) )
    U = as.vector( t(xx *eta)%*%(Phi  - tau) )
    U/cluster
  }
  
  
  DREfunc=function(L,R,T,x,delta,tau,ww,wr,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8);
    wl=(ww+wr)/(ww*wr); wl[is.nan(wl)]=0
    n=length(Y); nrl=sum(ifelse(delta==0,1,0));
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss =  pmax(1e-3, sqrt(diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wwind = ww*ind
    Phi = as.vector( pnorm( -res/ss ) * ww )
    U = as.vector( t(xx *eta)%*%(wwind - tau) )
    UR=matrix(0,p,1)
    UL=matrix(0,p,1)
    
    for (i in 1:n) {
      if(delta[i]==0){
        yindr=Y>=Y[i]
        denomr=sum(yindr*eta )*nrl
        Rft = as.vector( t(xx * eta *wr)%*%(as.numeric(( (yindr*eta*ind)-tau) *ind)) )
        UR=UR+(Rft/denomr)
      }
      if(delta[i]==0){
        yindl=Y<=Y[i]
        denoml=(sum((1-yindl) * eta )+1)*nrl
        Lft = as.vector( t(xx * eta*wl )%*%(as.numeric(( ((1-yindl)*eta*ind)-tau) *ind )) )
        UL=UL+((Lft/denoml))
      }
    }
    (U+(UR/nrl)+(UL/nrl))/(cluster)
  }
  
  
  Afunc=function(L,R,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y);
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    exactnum=sum(ifelse(delta==1,1,0))
    ss =  sqrt(pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    wwss=ww/ss
    phi = as.vector( dnorm( -res/ss )* (wwss))
    if(exactnum< (0.2*n)){phi = pmax(phi,0.05)}
    A = t(phi * xx  *eta) %*% xx + diag(p)*0.05
    A/cluster
  }
  
  Gfunc=function(L,R,x,delta,tau,ww,eta,cluster,beta,Sigma){
    L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y)
    xx=as.matrix(cbind(1,x)); p=ncol(xx)
    ss =  sqrt(pmax(1e-3, diag(xx%*%Sigma%*%t(xx))) ) 
    res = as.numeric(Y - xx%*%beta)
    ind = ifelse(res<=0,1,0)
    wwind = ww*ind
    
    Gam=xx*as.numeric(wwind-tau)
    Gamma=( t(Gam)%*%(Gam*eta) )
    GammaR=matrix(0,p,p)
    GammaL=matrix(0,p,p)
    
    for(i in 1:n){
      if(delta[i]==0){
        yindr=Y>=Y[i]
        denomr=sum(yindr*eta)
        nom=as.vector( t(xx*eta)%*% (yindr*wwind*eta) )
        R = nom/denomr
        GammaR=GammaR+(R)%*%t(R)
        
        yindl=Y<=Y[i]
        denoml=sum((1-yindl) * eta)+1
        nom=as.vector( t(xx*eta)%*% ((1-yindl)*wwind*eta) )
        L = nom/denoml
        GammaL=GammaL+(L)%*%t(L)
      }
    }
    (Gamma-GammaR+GammaL)/cluster
  }
  
  # # # update 'Gamma' matrix in sandwich variance (Zeng and Lin, 2008)
  Gfunc2 = function(L,R,x,delta,tau,ww,eta,cluster,beta,Sigma,B=100) {
    n=length(L); cennum=sum(ifelse(delta==1,1,0))
    library(MASS)
    Shat = t(replicate(B,{
      id = sample(n,n,replace = TRUE)
      Efunc(L=L[id],R=R[id],x=x[id,],delta=delta,tau=tau,ww=ww[id],eta=eta[id],cluster=cluster,beta = beta, Sigma = Sigma)
    }))
    Sigma = cov(Shat) * (cluster)
    Sigma
  }
  
  Gfunc3= function(L,R,x,delta,tau,ww,eta,id,cluster,beta,Sigma,B=100) {
    n=length(L)
    library(MASS)
    Shat = t(replicate(B,{
      tabid=as.vector(table(id))
      idx = as.vector(unlist(lapply(tabid, function(x) sample(x=x,size=x,replace = TRUE))))
      Efunc(L=L[idx],R=R[idx],x=x[idx,],delta=delta,tau=tau,ww=ww[idx],eta=eta[idx],cluster=cluster,beta = beta, Sigma = Sigma)
    }))
    Sigma = cov(Shat) * (cluster)
    Sigma
  }
  
  # update variance estimator
  up_Sigma = function(Y,Afunc,Gfunc,cluster){
    n=length(Y)
    invA = solve(Afunc)
    newSigma = (( (invA) %*% Gfunc %*% (invA) ) )
    newSigma/cluster
  }
  
  L = pmax(L,1e-8); Y=pmax(ifelse(delta==0,R,L),1e-8); n=length(Y); 
  if(is.null(id)){eta=rep(1,n); cluster=n}
  else{ci=rep(c(table(id)),c(table(id))); wi=(1/ci); eta=(wi^(index)); cluster=length(table(id))}
  
  if(wttype=="param"){ww=wtpicft(L=L,R=R,delta=delta);}
  else if(wttype=="Ishwaran" & n==sum(delta==1)){print("Use another weight estimating method (wttype=param or wttype=Beran)."); ww=NULL}
  else if(wttype=="Ishwaran"){ww=Ishrfwtpicfunc(L=L,R=R,delta=delta,x=x);}
  else if(wttype=="Beran" & is.null(hlimit)==F){ww=Berwtpicfunc(L=L,R=R,delta=delta,x=x,h=hlimit);}
  xx = as.matrix(cbind(1,x)); p = ncol(xx)
  old_beta = init = beta = PICrq(L=L,R=R,delta=delta,x=x,ww=ww,eta=eta,tau=tau)
  old_Sigma = Sigma = (diag(p)/cluster)
  
  i=0; eps=1; max.iter=100; tol = 1e-3; 
  while (i<max.iter & eps >= tol ) {
    Amat = Afunc(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    if(is.null(estimation)){
      new_beta = c(old_beta) - solve(Amat)%*%Efunc(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/(n)
    }
    else if(estimation=="dr"){
      wr=Rwtpicft(L=L,R=R,delta=delta)
      new_beta = c(old_beta) - solve(Amat)%*%DREfunc(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,wr=wr,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)/(n)
    }
    if(var.estimation=="IS" & is.null(id)==F){
      print("Use another variance estimating method (var.estimation=bootstrap)."); Gamma=NULL
    }
    else if(var.estimation=="IS"){
      Gamma = Gfunc(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }
    else if(var.estimation=="bootstrap" & is.null(id)){
      Gamma = Gfunc2(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,eta=eta,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
    }
    else if(var.estimation=="bootstrap"){
      Gamma = Gfunc3(L=L,R=R,x=x,delta=delta,tau=tau,ww=ww,eta=eta,id=id,cluster=cluster,beta = old_beta, Sigma = old_Sigma)
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
