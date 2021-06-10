library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)

set.seed(124)
### GENERATE DATA 
no <- 50
nx <- 20
dat <- coxed::sim.survdata(N = no, xvars = nx)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

# X <- matrix(rnorm(no * nx), nrow = no)
# yt <- rexp(no)
# yt[1] <- yt[2]
# #ys <- sample(0:1, size = no, replace = T)
# ys <- rep(1, no)
# 
# Xbig <- as.big.matrix(X)
# y <- cbind("time" = yt, "status" = ys)

##### CHECKLIST: CASES TO TEST #####
# no censored data, no ties (pass)
# censored data, no ties (pass)
# no censored data, ties (FAIL)
# censored data, ties (FAIL)
###################################

### MODEL PARAMETERS
alpha <- 0.5
penalty <- "enet" 
lambda  <- exp(seq(-1, -6, length.out = 8)) #c(0.1, 0.01)

##### FIT MODELS
fit.cx <- coxph(Surv(y) ~ X, ties = "breslow")
fit.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
fit.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)

beta.cx <- coef(fit.cx)
beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- as.matrix(fit.gn$beta)

beta <- beta.gn # used later when I try to reconstruct the deviance/likelihood values  

##### likelihoods
d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)

#satdev <- 2 * sum(d * log(d))
ll.sat <- -sum(d * log(d))
ll.cx  <- fit.cx$loglik

dev.cx <- -2 * (ll.cx - ll.sat)
dev.bl <- fit.bl$loss
dev.gn <- (1 - fit.gn$dev.ratio) * fit.gn$nulldev
dev.gn.null <- fit.gn$nulldev # defined as 2 * (ll_sat - ll_null)

# ll.gn.null <- ll.sat - dev.gn.null/2
# ll.gn.null
# dev.cx
# fit.gn$nulldev

calc.devs <- function(X, y, beta) {
  tOrder <- order(y[,1]) # indices of sorted times
  d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)
  dtime <- sort(unique(y[y[,2]==1,1]))  # sorted unique failure times for all uncensored obs.
  row.idx.cox <- which(y[tOrder,1] >= min(dtime)) # which (sorted) times >= the earliest failure time (should be all?)
  # NOTE: The above row.idx.cox omits the observations if they are censored AND are less than the earliest uncensored time
  # i.e. !(obs. are censored AND are less than the earliest uncensored time)
  d_idx <- integer(length(row.idx.cox)) # create empty vector of 0's of length row.idx.cox
  # index of the latest (unique) time <= the time of the i-th sorted obs.
  # in C++ code: "Index of unique failure time for subjects with failure; Index of the last unique failure time if censored"
  # QUESTION: Does the 'last' mean 'previous' in this context?
  for(i in 1:length(row.idx.cox)) d_idx[i] <- max(which(dtime <= y[tOrder[row.idx.cox[i]],1])) 
  
  tOrig <- order(tOrder) # map from the sorted times back to the original time ordering
  d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx) # unbelievably hacky solution to omitting any censored times occuring before the first uncensored one
  
  ### LOGLIKELIHOODS & DEVIANCES
  ll.sat <- -sum(d * log(d))
  
  ll <- 0
  for (j in 1:length(d)) {
    Dj <- d_idx2 == j # indices of the death/failures corresponding to the j-th unique failure time 
    Dj <- Dj[tOrig] & y[,2] # hacky solution to remove censored obs & place the indices in the original ordering
    term1 <- colSums(X[Dj,,drop=F]) %*% beta
    Rj <- (y[,1] >= y[Dj,1][1]) # all times in Dj will be the same, so we pick the first
    #term2 <- d[j] * log(sum(exp(X[Rj,,drop=F] %*% beta)))
    XRbeta <- X[Rj,,drop=F] %*% beta
    term2 <- d[j] * log(colSums(exp(XRbeta)))
    ll <- ll + (term1 - term2)
  }
  D <- -2 * (ll - ll.sat)
  
  ll.null <- 0
  for (j in 1:nrow(y)) {
    if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2] (no ties)
      Rj <- (y[,1] >= y[j,1])
      term2 <- log(sum(Rj))
      ll.null <- ll.null - term2
    }
  }
  D0 <- -2 * (ll.null - ll.sat)

  return (list("nulldev" = D0, "dev" = D, "ll_sat" = ll.sat, "ll" = ll, "ll_null" = ll.null))
}

# Dmat <- matrix(c("cx" = dev.cx, "bl" = c(NA, dev.bl), "gn" = c(dev.gn.null, dev.gn)), ncol = 2, nrow = 3, byrow = T)
# colnames(Dmat) <- c("D0", "D")
# rownames(Dmat) <- c("cx", "bl", "gn")

idx10 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] != 1 & y[,2] == 1)[1])
# two cens (diff times, one same as the uncens), one uncens
idx11 <- c(which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1])


idx <- idx11


pt <- proc.time()
devs1 <- calc.devs(X[idx,], y[idx,], beta.bl)
proc.time() - pt

coxnet.deviance(x = X[idx,], y = y[idx,], beta = beta.bl)

biglasso::cox.deviance(X = Xbig, y = y, beta = fit.gn$beta, row.idx = idx)$loglik
devs1$ll
