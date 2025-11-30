# Bound-constrained Nonstationary Gaussian Process
Tian Bai

November 30, 2025

This GitHub repository is the source code of the papaer **Bound-constrained Nonstationary
Gaussian Process
Regression for Ventilated Cavitation Prediction** by Tian Bai, Kuangqi Chen, and Dianpeng Wang, which is submitted to RESS.
This code allows the estimation and prediction for the Bound-constrained Nonstationary Gaussian Process (nbGP) model.

The paper and R package will be available soon.

The main function are built from R package `GPcluster` by Chih-Li Sung, and the details could be found in
https://github.com/ChihLi/GPcluster.

**The authors greatly thanks Sung et al. for their works.**

## Requirements
The nbGP model needs the following packages.
```r
library(GPcluster)
library(truncnorm)
library(SLHD) #For initial design of Section 4.2
library(tgp)
```

The main requirement is `GPcluster` package.
Please download the `GPcluster-main` folder in **THIS GitHub** and use the following code to install locally:
```r
devtools::install(".../GPcluster-master")
```
**PLEASE DO NOT INSTALL FROM THE ORIGINAL GPCluster GITHUB**

## Functions of nbGP-main.R
The main code of nbGP model is `nbGP-main.R`, including four main functions.
### Stratified sampling strategy from the truncated mixture distribution
```r
sample_truncated_mixture <- function(n, mu, sigma, l = -Inf, u = Inf)
```
- n: sample size;
- mu: predictive mean of $\mathbf{x}$;
- sigma: predictive standard variation of $\mathbf{x}$;
- l: lower bound at $\mathbf(x)$;
- u: upper bound at $\mathbf{x}$;
- Return: samples, the $1\times n$ samples from the truncated distribution for $\mathbf{x}$.

### Adjust the number of stratified sampling strategy in case the total number is not $n$
```r
adjust_matrix <- function(mat, target = 1000)
```
- mat: the matrix with each row defines the number of samples for each component in stratified
sampling strategy;
- target: the target number of samples;
- Return: adjusted_mat, the adjusted matrix with sum of each row equals target.

### Compute the truncated mean and variance using closed-form formulas
```r
calc_truncated_gaussian <- function(f, v, l, u)
```
- f: predictive mean of $\mathbf{x}$;
- v: predictive standard variation of $\mathbf{x}$;
- l: lower bound at $\mathbf{x}$;
- u: upper bound at $\mathbf{x}$;
- Return: list(ug = ug, vg2 = vg2), the truncated mean and variance of $\mathbf{x}$


### The main function of nbGP prediction.
```r
nbGP<-function(x_train,y_train,x_new,l_x,u_x,if_CI=FALSE,tilde_M=5,n_samples=1000,alpha=0.05,iter_max=100)
```
- x_train: training data $\mathbf{X}$;
- y_train: observations of tranining data $\mathbf{Y}$;
- x_new: new locations need to be predicted;
- l_x: the point-wise lower bounds of x_new;
- u_x: the point-wise upper bounds of x_new;
- if_Cl: default is `FALSE`, predict x_new using closed-form formulas; if `TRUE`, show the CI prediction using stratified sampling strategy;
- tilde_M: default is `5`, the maximum number of components;
- n_samples: default is `1000`; if `if_Cl==TRUE`, the sample size of stratified sampling strategy;
- alpha: default is `0.05`; if `if_Cl==TRUE`, return $1-\alpha$ Cl;
- iter_max: default is `100`; the maximum iteration for SEM algorithm.

## Example usage of Section 4.2
```r
testfun <- function(x) {
 sum1<-x[1]+1/4*x[2]
 if(x[1]>0.5&x[2]>0.5){
   return(exp(sum1))
 }else{
   return(1)}
}
lxfun<-function(x){
  lx<-ifelse(x[1]>0.5&x[2]>0.5,exp(0.5+1/4*0.5),0.8)
  return(lx)
}

uxfun<-function(x){
  ux<-ifelse(x[1]>0.5&x[2]>0.5,exp(x[1]+x[2]),1.2)
  return(ux)
}
n_samples <-1000
tilde_M<-5
alpha<-0.05
sdtrue<-0.01
n <-60
set.seed(as.numeric(Sys.time()))
D1<-maximinSLHD(t = 1, m = n,k =2)
x_train<-D1$StandDesign
y_train<-apply(x_train,1,testfun)
y_train <- y_train + rnorm(length(y_train),0, sd = sdtrue)
x1 <- runif(1000, min = 0,max= 1)
y1 <- runif(1000, min = 0, max=1)
x_new <- cbind(x1, y1)
y_new<-apply(x_new,1,testfun)
l_x<-apply(x_new, 1, lxfun)
u_x<-apply(x_new, 1, uxfun)
nbGP.fit<-nbGP(x_train,y_train,x_new,l_x,u_x,if_CI=FALSE,tilde_M,n_samples)
mu<-nbGP.fit$mu
var<-nbGP.fit$var
M_star<-nbGP.fit$M_star
```

## Citation
If you find our work helpful, feel free to give us a cite.
```
@misc{nbGP,
    title  = {Bound-constrained Nonstationary Gaussian Process Regression for Ventilated Cavitation Prediction},
    url    = {https://github.com/tbai114/Nonstationary-bounded-Gaussian-process},
    author = {Tian Bai, Kuangqi Chen, Dianpeng Wang},
    year   = {2025}
}
```

