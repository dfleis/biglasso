library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)
library(Rcpp)

set.seed(124)
### GENERATE DATA 
no <- 100
nx <- 50
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
lambda  <- exp(seq(-1, -6, length.out = 50)) #c(0.1, 0.01)

##### FIT MODELS
#fit.cx <- coxph(Surv(y) ~ X, ties = "breslow")
fit.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
fit.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)

beta.cx <- coef(fit.cx)
beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- as.matrix(fit.gn$beta)

##################### define deviance function fully in R for sanity tests
getDevsR <- function(X, y, beta, row.idx) {
  X <- X[row.idx,]
  y <- y[row.idx,]
  
  d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)
  ll.sat <- -sum(d * log(d))
  
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
  
  
  D_dR_sets <- lapply(1:length(d), function(j) {
    Dj <- d_idx2 == j
    Dj <- Dj[tOrig] & y[,2]
    dRj <- y[,1] == y[Dj, 1][1]
    list("Dj" = as.integer(which(Dj)),
         "dRj" = as.integer(which(dRj)))
  })
  
  expXRbeta_colsum <- 0 # rep(0, length(lambda)) # also, rep(0, ncol(beta)) works
  ll <- 0
  for (j in length(d):1) {
    Dj <- D_dR_sets[[j]]$Dj
    dRj <- D_dR_sets[[j]]$dRj
    term1 <- colSums(X[Dj,,drop=F]) %*% beta
    expXRbeta_colsum <- expXRbeta_colsum + colSums(exp(X[dRj,,drop=F] %*% beta))
    term2 <- d[j] * log(expXRbeta_colsum)
    ll <- ll + (term1 - term2)
  }
  ll <- as.numeric(ll)
  
  D <- -2 * (ll - ll.sat)
  return (list("D" = D, "ll" = ll, "ll.sat" = ll.sat))
}

##################### run calculations

row.idx <- seq(1, 100, by = 2)


system.time(out.r0 <- getDevsR(X = X, y = y, beta = fit.bl$beta, row.idx = 1:nrow(y)))
system.time(out.cpp0 <- cox.deviance(X = Xbig, y = y, beta = fit.bl$beta, row.idx = 1:nrow(y)))
system.time(out.gn0 <- glmnet::coxnet.deviance(y = y, x = X, beta = fit.bl$beta))

system.time(out.r1 <- getDevsR(X = X, y = y, beta = fit.bl$beta, row.idx = row.idx))
system.time(out.cpp1 <- cox.deviance(X = Xbig, y = y, beta = fit.bl$beta, row.idx = row.idx))
system.time(out.gn1 <- glmnet::coxnet.deviance(y = y[row.idx,], x = X[row.idx,], beta = fit.bl$beta))


################# figures (comparison)
plot(out.gn0 ~ lambda, log = 'x', type = 'l', lwd = 4, ylim = range(out.r0$D, out.cpp0$D, out.gn0))
lines(out.r0$D ~ lambda, col = 'red3', lwd = 4, lty = 'dashed')
lines(out.cpp0$D ~ lambda, col = 'deepskyblue3', lwd = 4, lty = 'dotted')

plot(out.gn1 ~ lambda, log = 'x', type = 'l', lwd = 4, ylim = range(out.r1$D, out.cpp1$D, out.gn1))
lines(out.r1$D ~ lambda, col = 'red3', lwd = 4, lty = 'dashed')
lines(out.cpp1$D ~ lambda, col = 'deepskyblue3', lwd = 4, lty = 'dotted')

plot(out.r1$D ~ lambda, log = 'x', type = 'l', col = 'red3', lwd = 4, ylim = range(out.r1$D, out.cpp1$D))
lines(out.cpp1$D ~ lambda, col = 'deepskyblue3', lwd = 4, lty = 'dotted')




out.r1$ll.sat
out.cpp1$ll.sat



