# Introduction to `ipcwqrPIC` R package


## Introduction
`ipcwqrPIC` is the R package to fit a new inverse-probability censoring weighted (IPCW) estimating procedure for censored quantile regression when the data are partially interval-censored that include doubly-censored (DC) data and partly interval-censored (PIC) and possibly correlated within same cluster.
Let $T$ and $X$ be the event time of interest and its related $p$-vector of covariates, respectively.
Our main objective is to estimate 
the $p$-dimensional quantile coefficient vector ${\boldsymbol{\beta}}_0(\tau)$
for some $\tau \in[\tau_L,\tau_R]\subset (0, 1)$ 
in the following linear quantile regression model:

$$T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0(\tau) + e_i(\tau),\quad i=1, \ldots ,n, $$

where $e_i(\tau)$ is the random error 
whose $\tau$ th quantile conditional on 
${\bf x}_i$ equals 0. 
When the data are subject to partially interval-censoring, 
left and right endpoints of the censoring time, $L$ and $R$,
are observed instead of $T$ such that $T\in(L,R)$.
Note that double-censoring  can also  be viewed as 
a special case of partly interval-censoring, 
i.e., $T$ is left-censored if $L=0$ and right-censored if $R=\infty$. 


## Description
This R package `ipcwqrPIC` implements an inverse-probability censoring weighted (IPCW) procedure for censored quantile regression for (cluster-correlated) partially interval-censored data, which includes both double-censoring and partially interval-censoring.

Vignettes is available in [here](http://htmlpreview.github.io/?https://github.com/YejiStat/ipcwqrPIC/blob/main/vignettes/ipcwqrPIC.html).


## Usages 
```{r message=FALSE, warning=FALSE}
library(PICBayes)
data("mCRC")
d = with(data.frame(mCRC), data.frame(U = ifelse(y==0,R,L),
                                      V = ifelse(y==2,L,R),
                                      # Cluster weighted data
                                      id=(rep(c(table(SITE)),c(table(SITE)))),
                                      # Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.
                                      x1= case_when(TRT_C == 0 ~ 0, #Pan et al data
                                                    TRT_C == 1 ~ 1),
                                      # Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.
                                      x2= case_when(KRAS_C == 0 ~ 1,
                                                    KRAS_C == 1 ~ 0),
                                      y = y
)); table(d$y)
delta=with(d,ifelse(y==3,1,
             ifelse(y==2,2,
                           ifelse(y==1,4,3)))); table(delta)
#
#   0   1   2   3 
# 168 329 306  52
delta=with(d,ifelse(y==3,1,
             ifelse(y==2,2,
                           ifelse(y==1,4,3)))); table(delta)
# delta
#   1   2   3   4 
#  52 306 168 329
U=(log(d$U));V=log(d$V); x = cbind(d$x1,d$x2); id=d$id;  tau=0.3;
ipcwqrPIC::picrq(L=U,R=V,delta=delta,x=x,tau=tau,wttype="KM",application=TRUE,var.estimation = "IS",id=id,index = 1)
#           tau coefficients       se   pvalue  lower bd upper bd
# Intercept 0.3     3.748561 0.535708 0.000000  2.698573 4.798549
# 2         0.3    -0.111274 0.488262 0.409863 -1.068267 0.845719
# 3         0.3     0.294631 0.547730 0.295318 -0.778919 1.368181
ipcwqrPIC::picrq(L=U,R=V,delta=delta,x=x,tau=tau,wttype="Beran",application=TRUE,hlimit=0.1,var.estimation = "Bootstrap",id=id,index = 1,B=100)
#           tau coefficients       se   pvalue  lower bd upper bd
# Intercept 0.3     2.707818 0.177515 0.000000  2.359888 3.055748
# 2         0.3     0.082049 0.215093 0.351432 -0.339534 0.503631
# 3         0.3     1.141851 0.222131 0.000000  0.706475 1.577226
ipcwqrPIC::picrq(L=U,R=V,delta=delta,x=x,tau=tau,wttype="Beran",application=TRUE,estimatio="DR",hlimit=0.1,var.estimation = "Bootstrap",id=id,index = 1,B=100)
#           tau coefficients       se   pvalue  lower bd upper bd
# Intercept 0.3     2.719725 0.176441 0.000000  2.373901 3.065549
# 2         0.3     0.081384 0.219507 0.355409 -0.348849 0.511617
# 3         0.3     1.166419 0.207143 0.000000  0.760419 1.572420
```


## References


* Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.

* Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.

* Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time-
dependent covariates. Biostatistics. 18 (1), 132–146.

* Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.

* Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). 
Inverse weighted quantile regression with partially interval-censored data.
*Submitted to BMJ*.
