library(MASS)
library(splines)
library(latex2exp)
source('../toSource.R')

expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

n <- 4*2000
nsim <- 500
rateseq <- seq(0.1,0.5,by=0.05)

set.seed(1)

# Dimensions to try
dims <- c(1,5,10,50)

for (D in c(4)) {
  for (d in dims) {
    res2 <- cate_hat <- NULL
    
    for (rate in rateseq) {
      print(paste("Rate:", rate, "Dim:", d))
      res <- data.frame(matrix(nrow=nsim, ncol=6))
      colnames(res) <- c("forster_poly", "ls_poly","forster_ns","ls_ns", "forster_bs", "ls_bs")
      cate_hat_mat <- data.frame(matrix(nrow=nsim, ncol=6))
      
      for (i in 1:nsim) {
        if (i %% 100 == 0) print(i)
        
        s <- sort(rep(1:4, n/4))
        
        # Multidimensional X
        x <- matrix(0.5*rnorm(n*d) + runif(n*d), nrow=n, ncol=d)
        
        # More complex DGP
        # propensity score depends on first two dimensions nonlinearly
        ps <- 0.1 + 0.8 * (x[,1] > 0)
        ps <- pmin(pmax(ps, 0.01), 0.99) # keep in (0,1)
        
        # potential outcomes: additive + interaction + nonlinear
        mu0 <- (x[,1] <= -.5)*0.5*(x[,1]+2)^2 + (x[,1]/2+0.875)*(x[,1]>-1/2 & x[,1]<0) +
          (x[,1]>0 & x[,1]<.5)*(-5*(x[,1]-0.2)^2 +1.075) + (x[,1]>.5)*(x[,1]+0.125);        
        mu1 <- mu0
        tau <- 0
        
        # treatment assignment
        a <- rbinom(n, 1, ps)
        y <- a*mu1 + (1-a)*mu0 + rnorm(n, sd=0.2 - 0.1*cos(2*pi*x[,1]))

        ## estimate nuisance functions
        # for simplicity, use only first covariate in smoothing
        pihat <- expit(logit(ps) + rnorm(n, mean=1/(n/4)^rate, sd=1/(n/4)^rate))
        mu1hat <- predict(smooth.spline(x[a==1 & s==2,1], y[a==1 & s==2], cv=T), x[,1])$y
        mu0hat <- predict(smooth.spline(x[a==0 & s==2,1], y[a==0 & s==2], cv=T), x[,1])$y
        
        ## construct estimators
        plugin <- mu1hat - mu0hat
        pseudo <- ((a - pihat) / (pihat*(1-pihat))) * (y - a*mu1hat - (1-a)*mu0hat) + mu1hat - mu0hat
        
        tau_forster_poly <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="forster", basis_type="poly")
        tau_forster_ns   <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="forster", basis_type="ns")
        tau_forster_bs   <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="forster", basis_type="bs")
        
        tau_ls_poly <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="ls", basis_type="poly")
        tau_ls_ns   <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="ls", basis_type="ns")
        tau_ls_bs   <- series_df(x[s==3,1], pseudo[s==3], x[s==4,1], df=D, type="ls", basis_type="bs")
        
        test_idx <- which(s==4)
        res$forster_poly[i] <- (n/4)*mean((tau_forster_poly[[1]] - tau)^2)
        res$forster_ns[i]   <- (n/4)*mean((tau_forster_ns[[1]] - tau)^2)
        res$forster_bs[i]   <- (n/4)*mean((tau_forster_bs[[1]] - tau)^2)
        res$ls_poly[i]      <- (n/4)*mean((tau_ls_poly[[1]] - tau)^2)
        res$ls_ns[i]        <- (n/4)*mean((tau_ls_ns[[1]] - tau)^2)
        res$ls_bs[i]        <- (n/4)*mean((tau_ls_bs[[1]] - tau)^2)
        
        cate_hat_mat$forster_poly[i] <- mean(tau_forster_poly[[1]])
        cate_hat_mat$forster_ns[i]   <- mean(tau_forster_ns[[1]])
        cate_hat_mat$ls_poly[i]      <- mean(tau_ls_poly[[1]])
        cate_hat_mat$ls_ns[i]        <- mean(tau_ls_ns[[1]])
      }
      
      res2 <- rbind(res2, c(rate, apply(res,2,mean)))
      cate_hat <- rbind(cate_hat, c(rate, apply(cate_hat_mat,2,mean)))
    }
    
    save(res2, file = paste0('mixture_forster_ls_', n/4, '_d', d, '_new.RData'))
  }
}
