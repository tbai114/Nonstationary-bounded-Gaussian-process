library(GPcluster)
library(truncnorm)
library(SLHD)
library(tgp)

sample_truncated_mixture <- function(n, mu, sigma, l = -Inf, u = Inf) {
  is_l_finite <- is.finite(l)
  is_u_finite <- is.finite(u)
  if (!is_l_finite) {
    Z <- pnorm((u - mu)/sigma)
  }
  else if (!is_u_finite) {
    Z<- 1-pnorm((l - mu)/sigma)
    }
  else {
      Z <- pnorm((u - mu)/sigma) - pnorm((l -mu)/sigma)
  }
  samples <- numeric(n)
  for (i in 1:n) {
    if (runif(1) < Z) {
      samples[i] <- rtruncnorm(1, a = l, b = u, mean = mu, sd = sigma)
    } else {
      if (!is_l_finite && !is_u_finite) {
        samples[i] <- NA_real_
      } else if (!is_l_finite) {
        samples[i] <- u
      } else if (!is_u_finite) {
        samples[i] <- l
      } else {
        prob_l <- pnorm((l - mu)/sigma) / (1 - Z)
        if (runif(1) < prob_l) {
          samples[i] <- l
        } else {
          samples[i] <- u
        }
      }
    }
  }
  return(samples)
}

adjust_matrix <- function(mat, target = 1000) {
  adjusted_mat <- mat
  for (i in 1:nrow(adjusted_mat)) {
    row <- adjusted_mat[i, ]
    current_sum <- sum(row)
    if (current_sum == target) {
      next
    }
    non_zero <- which(row != 0)
    if (length(non_zero) == 0) {
      warning(paste("Row", i, "is all zero, skipping adjustment."))
      next
    }
    max_non_zero_col <- max(non_zero)
    delta <- target - current_sum
    adjusted_mat[i, max_non_zero_col] <- adjusted_mat[i, max_non_zero_col] + delta
  }
  return(adjusted_mat)
}

calc_truncated_gaussian <- function(f, v,l, u){
  is_l_finite <- is.finite(l)
  is_u_finite<-is.finite(u)
  if (is_l_finite && is_u_finite) {
    print('try')
    Z <- pnorm(u, mean = f, sd = v) - pnorm(l, mean = f, sd = v)
    alpha <- (l - f) / v
    beta <- (u - f) / v
    ug <- Z * f + (dnorm(alpha) - dnorm(beta)) * v + l * pnorm(alpha) + u * (1 - pnorm(beta))
    print(ug)
    vg2 <- Z * (v^2 + f^2) + 2 * v * f * (dnorm(alpha) - dnorm(beta)) + v^2 * (alpha * dnorm(alpha) - beta * dnorm(beta)) + 
      l^2 * pnorm(alpha) + 
      u^2 * (1 - pnorm(beta)) - 
      ug^2
  }
  else if (!is_u_finite){
        Z<-1-pnorm(l, mean = f, sd = v)
        alpha <- (l-f)/v
        ug <- Z * f + dnorm(alpha) * v+l* pnorm(alpha)
        vg2 <- Z* (v^2 + f^2) +
          2 * v * f * dnorm(alpha)+
          v^2 * alpha * dnorm(alpha) +
          l^2* pnorm(alpha) -ug^2
        }
  else if (!is_l_finite) {
          Z <- pnorm(u, mean = f, sd = v)
          beta <- (u - f) /v
          ug <- Z *f - dnorm(beta) * v + u * (1- pnorm(beta))
          vg2 <- Z * (v^2 + f^2)-
          2 * v * f * dnorm(beta) -
            v^2 * beta * dnorm(beta) +
            u^2 * (1 - pnorm(beta))-
            ug^2
  }
  return(list(ug = ug, vg2 = vg2))
}

nbGP<-function(x_train,y_train,x_new,l_x,u_x,if_CI=FALSE,tilde_M=5,n_samples=1000,alpha=0.05,iter_max=100){
   n<-length(x_train[,1])
   ############Fit the nugget using tgp##########
   train_idx<- sample(1:n, size = floor(0.7*n))
   test_idx <- setdiff(1:n, train_idx)
   X_train <-x_train[train_idx,]
   Z_train <- y_train[train_idx]
   X_test <-x_train[test_idx,]
   Z_test<-y_train[test_idx]
   fit<-btgp(X_train,Z_train,X_test)
   Z_pre<-fit$ZZ.km
   Z_var<-fit$ZZ.ks2
   nugget<-mean((Z_test-Z_pre)^2)-mean(Z_var)
   sd_train<-sd(y_train)
   nugget_1<-nugget/(sd_train^2)
   ############Fit the nbGP model################
   LOOCV_min<-1e6
   fit.object<-0
   M<-2
   M_star<-1
   while(M<=tilde_M){
     try(fit.object_try <-GPcluster_fit(x_train, y_train, K = M, iter_max = iter_max,nugget=nugget_1),
         silent = TRUE
     )
     if (inherits(fit.object_try, "try-error")) {
       next
     }
     LOOCV_M<-min(fit.object_try$LOOCV)
     if(LOOCV_M<LOOCV_min){
       fit.object<-fit.object_try
       M_star<-M
       LOOCV_min<-LOOCV_M
     }
     M<-M+1
   }
   pred.out <- predict(fit.object, x_new,conf.out=1-alpha)
   yhat<-pred.out$y_ori
   yhat<-yhat*c(sd(y_train)) + mean(y_train)
   sehat<-pred.out$se_ori
   sehat<-sqrt(sehat^2*c(var(y_train)))
   prob<-pred.out$prob
   if(if_CI==TRUE){
     #Compute the mean, variance, and Cl using stratified sampling#
     sample_counts <- round(prob * n_samples)
     sample_counts<-adjust_matrix(sample_counts,n_samples)
     samples_matrix <- matrix(0, nrow = length(x_new[,1]), ncol = n_samples)
     for (i in 1:length(x_new[,1])) {
       l=l_x[i]
       u = u_x[i]
       all_samples <- numeric(0)
       for(k in 1:M_star){
         mu=yhat[i,k]
         sigma=sehat[i,k]
         n_samples_k<-sample_counts[i,k]
         n_samples_final<-sample_counts[i,M_star]
         if (n_samples_k > 0){
           samples<-sample_truncated_mixture(n_samples_k, mu, sigma, l, u)
           all_samples <- c(all_samples, samples)
         }
       }
         samples_matrix[i, ]<-all_samples
       }
   mu<- apply(samples_matrix, 1, mean)
   var<- apply(samples_matrix, 1, var)
   ci_lower<- apply(samples_matrix, MARGIN = 1, function(x) quantile(x, probs = alpha/2))
   ci_upper<-apply(samples_matrix, MARGIN = 1, function(x) quantile(x, probs = 1-alpha/2))
   return(list(mu=mu,var=var,M_star=M_star,ci_lower=ci_lower,ci_upper=ci_upper))
   }
   else{
     #Compute the mean and variance with closed-form formulas
     uall<-matrix(NA,length(x_new[,1]),ncol=M_star, nrow=length(x_new[,1]))
     vall<-matrix(NA,length(x_new[, 1]),ncol=M_star,nrow=length(x_new[,1]))
     for (i in 1:length(x_new[,1])) {
       l=l_x[i]
       u=u_x[i]
       for(k in 1:M_star){
        mu_x=yhat[i,k]
        sigma_x=sehat[i,k]
        re1<-calc_truncated_gaussian(mu_x,sigma_x,l,u)
        uall[i,k]<-re1$ug
        vall[i,k]<-re1$vg2
       }
     }
    mu<-apply(prob*uall, 1,sum)
    var<-apply(prob * (uall^2+vall), 1,sum)-mu^2
    return(list(mu=mu,var=var,M_star=M_star))
   }
}

####Example usage of Section 4.2
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

