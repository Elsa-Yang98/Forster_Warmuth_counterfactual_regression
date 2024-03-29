rm(list = ls())
library(MASS)
library(splines)
library(latex2exp)
source('toSource.R')

D = 4

expit <- function(x){ exp(x)/(1+exp(x)) };

logit <- function(x){ log(x/(1-x)) }

n <- 4*2000; nsim <- 500;

rateseq <- seq(0.1,0.5,by=0.05);


cate_hat <- res2 <- NULL
set.seed(1)
for (rate in rateseq){
  
  print(rate)
  
  cate_hat_mat <- data.frame(matrix(nrow=nsim,ncol=14))
  colnames(cate_hat_mat) <- c("plugin","xl","drl","drl_poly","oracle.drl","bs","ns","poly","forster_poly","forster_ns", "ls_poly","ls_ns", "forster_bs", "ls_bs")
  res <- data.frame(matrix(nrow=nsim,ncol=14))
  colnames(res) <- c("plugin","xl","drl","drl_poly","oracle.drl","bs","ns","poly","forster_poly","forster_ns", "ls_poly","ls_ns", "forster_bs", "ls_bs")
  for (i in 1:nsim){
    if(i %%100==0){print(i)}
    s <- sort(rep(1:4,n/4));
    x <- (runif(n,-1,1));
    ps <- 0.1 + 0.8*(x>0)
    mu0 <- (x <= -.5)*0.5*(x+2)^2 + (x/2+0.875)*(x>-1/2 & x<0) +
      (x>0 & x<.5)*(-5*(x-0.2)^2 +1.075) + (x>.5)*(x+0.125);
    mu1 <- mu0;
    tau <- 0
    a <- rbinom(n,1,ps);
    y <- a*mu1 + (1-a)*mu0 + rnorm(n,sd=(.2-.1*cos(2*pi*x)))
    
    ## estimate nuisance functions
    
    pihat <- expit( logit(ps) + rnorm(n,mean=1/(n/4)^rate,sd=1/(n/4)^rate))
    
    mu1hat <- predict(smooth.spline(x[a==1 & s==2],y[a==1 & s==2]),x)$y
    mu0hat <- predict(smooth.spline(x[a==0 & s==2],y[a==0 & s==2]),x)$y
    
    ## construct estimators
    
    plugin <- mu1hat-mu0hat
    
    x1 <- predict(smooth.spline(x[a==1 & s==3],(y-mu0hat)[a==1 & s==3]),x)$y
    x0 <- predict(smooth.spline(x[a==0 & s==3],(mu1hat-y)[a==0 & s==3]),x)$y
    xl <- pihat*x0 + (1-pihat)*x1
    
    pseudo <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
    drl <- predict(smooth.spline(x[s==3],pseudo[s==3]),x)$y
    
  
    
    tau_forster_poly = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "poly")
    tau_forster_ns = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "ns")
    tau_forster_bs = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "bs")
    
    tau_ls_poly = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "poly")
    tau_ls_ns = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "ns")
    tau_ls_bs = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "bs")
    
    
    pseudo.or <- ((a-ps)/(ps*(1-ps)))*(y-a*mu1-(1-a)*mu0) + mu1-mu0
    oracle.drl <- predict(smooth.spline(x[s==3],pseudo.or[s==3]),x)$y
    ## save MSEs
    res$plugin[i] <- (n/4)*mean((plugin-tau)[s==4]^2)
    res$xl[i] <- (n/4)*mean((xl-tau)[s==4]^2)
    res$drl[i] <- (n/4)*mean((drl-tau)[s==4]^2)
    #res$drl_poly[i] <- (n/4)*mean((drl_hat_poly-tau)[s==4]^2)       
    #res$drl_ns[i] <- (n/4)*mean((drl_hat_ns-tau)[s==4]^2)
    res$oracle.drl[i] <- (n/4)*mean((oracle.drl-tau)[s==4]^2)
    res$forster_poly[i] <- (n/4)*mean((tau_forster_poly[[1]]-tau)^2)
    res$forster_ns[i] <- (n/4)*mean((tau_forster_ns[[1]]-tau)^2)
    res$forster_bs[i] <- (n/4)*mean((tau_forster_bs[[1]]-tau)^2)
    res$ls_poly[i] <-(n/4)*mean((tau_ls_poly[[1]]-tau)^2)
    res$ls_ns[i] <- (n/4)*mean((tau_ls_ns[[1]]-tau)^2)
    res$ls_bs[i] <- (n/4)*mean((tau_ls_bs[[1]] -tau)^2)
    
    cate_hat_mat$plugin[i] <- mean(plugin)
    cate_hat_mat$xl[i] <- mean(xl)
    cate_hat_mat$drl[i] <- mean(drl)
    #cate_hat_mat$drl_poly[i] <- mean(drl_hat_poly)
    #cate_hat_mat$drl_ns[i] <- mean(drl_hat_ns)
    cate_hat_mat$oracle.drl[i] <- mean(oracle.drl)
    cate_hat_mat$forster_poly[i] <- mean(tau_forster_poly[[1]])
    cate_hat_mat$forster_ns[i] <- mean(tau_forster_ns[[1]])
    cate_hat_mat$forster_bs[i] <- mean(tau_forster_bs[[1]] )
    cate_hat_mat$ls_poly[i] <-mean(tau_ls_poly[[1]])
    cate_hat_mat$ls_ns[i] <- mean(tau_ls_ns[[1]])
    cate_hat_mat$ls_bs[i] <- mean(tau_ls_bs[[1]] )
    
  }
  res2 <- rbind(res2, c(rate, apply(res,2,mean)))
  cate_hat<- rbind(cate_hat, c(rate, apply(cate_hat_mat,2,mean)))
}

final = list(res2, cate_hat)


save(final, file = "plotdata_uniform_2000.RData")










