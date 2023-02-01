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
```{r}
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
                                      delta = case_when(IC == 0 ~ 1,
                                                        IC == 1 ~ 0)
));
L=(log(d$U));R=log(d$V); delta=d$delta
x = cbind(d$x1,d$x2); id=d$id;  tau=0.1;
ipcwqrPIC::picrq(L,R,delta,x=x,tau=tau,var.estimation = "IS")
#>           tau coefficients       se   pvalue  lower bd upper bd
#> Intercept 0.1     2.548123 0.115091 0.000000  2.322544 2.773702
#> 2         0.1     0.170877 0.151378 0.129488 -0.125823 0.467577
#> 3         0.1     0.240066 0.149144 0.053739 -0.052256 0.532387
ipcwqrPIC::picrq(L,R,delta,x=x,tau=tau,var.estimation = "IS",id=id,index = 1)
#>           tau coefficients       se   pvalue  lower bd upper bd
#> Intercept 0.1     2.549014 0.096420 0.000000  2.360030 2.737998
#> 2         0.1     0.167647 0.125468 0.090746 -0.078270 0.413563
#> 3         0.1     0.248939 0.123488 0.021905  0.006903 0.490975
ipcwqrPIC::picrq(L,R,delta,x=x,tau=tau,wttype="Beran",hlimit=0.9,var.estimation = "IS",id=id,index = 1)
#>           tau coefficients       se   pvalue  lower bd upper bd
#> Intercept 0.1     2.548473 0.096867 0.000000  2.358614 2.738332
#> 2         0.1     0.167874 0.126269 0.091841 -0.079612 0.415361
#> 3         0.1     0.248816 0.124234 0.022600  0.005317 0.492315
ipcwqrPIC::picrq(L,R,delta,x=x,tau=tau,estimation = "DR",var.estimation = "IS",id=id,index = 1)
#>           tau coefficients       se   pvalue  lower bd upper bd
#> Intercept 0.1     2.400826 0.097163 0.000000  2.210387 2.591265
#> 2         0.1     0.134126 0.116797 0.125409 -0.094797 0.363049
#> 3         0.1     0.240651 0.115571 0.018658  0.014133 0.467169
```


## References


* Beran, R. (1981). Nonparametric Regression with Randomly Censored Survival Data. Technical Report, Univ.California, Berkeley.

* Ishwaran, H., U. B. Kogalur, E. H. Blackstone, and M. S. Lauer (2008). Random survival forests. The annals of applied statistics. 2 (3), 841–860.

* Chiou, S. H., Kang, S. and Yan, J. (2015). 
Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach.
Statistics in Medicine 34(9): 1495–-1510.

* Gorfine, M., Y. Goldberg, and Y. Ritov (2017). A quantile regression model for failure-time data with time-
dependent covariates. Biostatistics. 18 (1), 132–146.

* Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating equations with general weight for accelerated failure time models: an induced smoothing approach. Statistics in Medicine 34(9): 1495–-1510.

* Zeng, D. and D. Lin (2008). Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9 (2), 355–363.

* Pan, C. (2021). 
PICBayes: Bayesian Models for Partly Interval-Censored Data. R package. 
https://CRAN.R-project.org/package=PICBayes.

* Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D. (2023+). 
Inverse weighted quantile regression with partially interval-censored data.
*Submitted to SMMR*.
