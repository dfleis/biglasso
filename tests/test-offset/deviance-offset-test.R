library(microbenchmark)
library(survival)
library(glmnet)

set.seed(124)
n <- 250
p <- 11
X <- matrix(rnorm(n * p), nrow = n)
y <- cbind("time" = round(rexp(n)*5)+1, "status" = rbinom(n, 1, 0.75))
offset <- rnorm(n)

fit.cx <- coxph(Surv(y) ~ X + offset(offset), ties = "breslow")
fit.gn <- glmnet(x = X, y = y, offset = offset, family = "cox", lambda = c(0, 0.01))

d      <- as.numeric(table(y[y[,2]==1,1]))
llsat <- -sum(d * log(d))



dev.gn     <- coxnet.deviance(x = X, y = y, offset = offset, beta = coef(fit.gn))
nulldev.gn <- coxnet.deviance(x = X, y = y, offset = offset, beta = coef(fit.gn)*0)
ll.gn     <- -dev.gn/2 + llsat
llnull.gn <- -nulldev.gn/2 + llsat
devs.beta.gn <- rbind("glmnet" = c("null" = llnull.gn, "fit" = ll.gn), "coxph" = fit.cx$loglik)

dev.gn     <- coxnet.deviance(x = X, y = y, offset = offset, beta = coef(fit.cx))
nulldev.gn <- coxnet.deviance(x = X, y = y, offset = offset, beta = coef(fit.cx)*0)
ll.gn     <- -dev.gn/2 + llsat
llnull.gn <- -nulldev.gn/2 + llsat
devs.beta.cx <- rbind("glmnet" = c("null" = llnull.gn, "fit" = ll.gn), "coxph" = fit.cx$loglik)

devs.beta.gn
devs.beta.cx

beta <- coef(fit.cx)
Xbig <- bigmemory::as.big.matrix(X)




########################################################################################
########################################################################################

loglik1 <- function(XX, yy, bb, offs) {
  ll <- 0
  for (i in 1:nrow(yy)) {
    C <- yy[i,2]
    term1 <- XX[i,] %*% bb + offs[i]
    
    Rset <- yy[,1] >= yy[i,1]
    term2 <- log(sum(exp(XX[Rset,,drop=F] %*% bb + offs[Rset])))
    
    ll <- ll + C * (term1 - term2)
  }
  return (ll)
}

loglik2 <- function(XX, yy, bb, offs) {
  ll <- 0
  for (i in 1:nrow(yy)) {
    if (yy[i,2] == 1) { 
      term1 <- XX[i,] %*% bb + offs[i]
      
      Rset <- yy[,1] >= yy[i,1]
      term2 <- log(sum(exp(XX[Rset,,drop=F] %*% bb + offs[Rset])))
    
      ll <- ll + term1 - term2
    }
  }
  return (ll)
}

loglik3 <- function(XX, yy, bb, offs) {
  d     <- as.numeric(table(yy[yy[,2]==1,1]))
  dtime <- sort(unique(yy[yy[,2]==1,1]))
  
  ll <- 0
  for (j in 1:length(d)) {
    Dset <- (yy[,1] == dtime[j]) & (yy[,2] == 1)
    S <- colSums(XX[Dset,,drop=F])
    term1 <- S %*% bb + sum(offs[Dset])
    
    Rset <- (yy[,1] >= dtime[j]) #| Dset
    term2 <- d[j] * log(sum(exp(XX[Rset,,drop=F] %*% bb + offs[Rset])))
    
    ll <- ll + term1 - term2
  }
  return (ll)
}

loglik.null1 <- function(yy, offs) {
  ll <- 0
  for (i in 1:nrow(yy)) {
    C <- yy[i,2]
    term1 <- offs[i]
    
    Rset <- yy[,1] >= yy[i,1]
    term2 <- log(sum(exp(offs[Rset])))
    
    ll <- ll + C * (term1 - term2)
  }
  return (ll)
}

loglik.null2 <- function(yy, offs) {
  ll.null <- 0
  for (j in 1:nrow(yy)) {
    if (yy[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2]
      R <- yy[,1] >= yy[j,1]
      ll.null <- ll.null + offs[j] - log(sum(exp(offs[R])))
    }
  }
  return (ll.null)
}

loglik.null3 <- function(yy, offs) {
  d     <- as.numeric(table(yy[yy[,2]==1,1]))
  dtime <- sort(unique(yy[yy[,2]==1,1]))
  
  ll.null <- 0
  for (j in 1:length(d)) {
    Dset <- (yy[,1] == dtime[j]) & (yy[,2] == 1)
    term1 <- sum(offs[Dset])
    
    Rset <- (yy[,1] >= dtime[j]) #| Dset
    term2 <- d[j] * log(sum(exp(offs[Rset])))
    
    ll.null <- ll.null + term1 - term2
  }
  return (ll.null)
}

loglik.tmp <- function(XX, yy, bb, offs) {
  d     <- as.numeric(table(yy[yy[,2]==1,1]))
  dtime <- sort(unique(yy[yy[,2]==1,1]))
  
  ### fitted likelihoods
  tOrder      <- order(yy[,1])
  row.idx.cox <- which(yy[tOrder,1] >= min(dtime))
  d_idx       <- integer(length(row.idx.cox))
  for(i in 1:length(row.idx.cox)) d_idx[i] <- max(which(dtime <= yy[tOrder[row.idx.cox[i]],1]))
  
  tOrig <- order(tOrder) # map from the sorted times back to the original time ordering
  # unbelievably hacky solution to omitting any censored times occurring before the first uncensored one:
  d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx)

  D_dR_sets <- lapply(1:length(d), function(j) {
    dRset <- (d_idx2 == j)[tOrig]
    Dset  <- dRset & y[,2]
    list("D"  = as.integer(which(Dset) - 1 + 1), 
         "dR" = as.integer(which(dRset) - 1 + 1))
  })
  
  ll <- 0
  expXRbeta_colsum <- rep(0, ncol(bb))
  for (j in length(D_dR_sets):1) {
    Dset  <- D_dR_sets[[j]]$D
    dRset <- D_dR_sets[[j]]$dR
    
    S     <- colSums(XX[Dset,,drop=F])
    term1 <- S %*% bb + sum(offs[Dset])
    
    expXRbeta_colsum <- expXRbeta_colsum + colSums(exp(XX[dRset,,drop=F] %*% bb + offs[dRset]))
    term2 <- d[j] * log(expXRbeta_colsum)
    
    ll <- ll + term1 - term2
    
    #print(sprintf("expXRbeta_colsum[%i] = %f", j, as.numeric(expXRbeta_colsum)))
    print(sprintf("term2[%i] = %f", j, as.numeric(term2)))
  }
  return (ll)
}

########################################################################################
########################################################################################

btmp <- coef(fit.cx)
b <- as(matrix(btmp), "dgTMatrix")

out1 <- loglik.tmp(XX = X, yy = y, bb = b, offs = offset)
out2 <- biglasso::cox.deviance(X = Xbig, y = y, beta = b, offset = offset)
out1
out2$loglike
out2$loglike
fit.cx$loglik
out1

out$loglike



out <- loglik.tmp(XX = X, yy = y, bb = b, offs = offset)
out

pt <- proc.time()
mb <- microbenchmark(loglik1(XX = X, yy = y, bb = beta, offs = offset),
                     loglik2(XX = X, yy = y, bb = beta, offs = offset),
                     loglik3(XX = X, yy = y, bb = beta, offs = offset),
                     loglik.null1(yy = y, offs = offset),
                     loglik.null2(yy = y, offs = offset),
                     loglik.null3(yy = y, offs = offset), times = 10, unit = "ms")
proc.time() - pt
boxplot(mb, outline = F)
summary(mb)

#glmnet::coxnet.deviance(x = X, y = y, offset = offset, beta = b*0)



# library(Rcpp)
# library(RcppArmadillo)
# 
# cppFunction(depends = "RcppArmadillo", '
#   void f(NumericVector v) {
#     arma::mat m = arma::mat(1, v.size());
#     for (int i = 0; i < m.size(); i++) {
#       m[i] = v[i];
#     }
#     
#     for (int i = 0; i < m.size(); i++) {
#       Rprintf("m[%i] = %f\\n", i, m[i]);
#     }
#   }
# ')
# 
# x <- rnorm(5)
# f(x)

















