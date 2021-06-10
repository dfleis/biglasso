library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)
library(Rcpp)

set.seed(124)
### GENERATE DATA 
no <- 10000
nx <- 100
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
# no censored data, ties 
# censored data, ties 
###################################

### MODEL PARAMETERS
alpha <- 0.5
penalty <- "enet" 
lambda  <- exp(seq(-1, -6, length.out = 100)) #c(0.1, 0.01)

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

###### manual tests 
ll.sat <- -sum(d * log(d))

# D_R_sets <- lapply(1:length(d), function(j) { 
#   # NOTE: structure is slightly different in the deviance.cox function
#   # so that it's more amenable to passing it to C++, i.e.
#   # indexes starting at 0, unnecessary inclusion of d[j].
#   Dj <- d_idx2 == j
#   Dj <- Dj[tOrig] & y[,2]
#   Rj <- (y[,1] >= y[Dj,1][1])
#   list("Dj" = as.integer(which(Dj)), 
#        "Rj" = as.integer(which(Rj)),
#        "dj" = d[j]) 
# })
# ll.list <- sapply(D_R_sets, function(DRs) {
#   term1 <- colSums(X[DRs$Dj,,drop=F]) %*% beta
#   term2 <- DRs$d * log(colSums(exp(X[DRs$Rj,,drop=F] %*% beta)))
#   term1 - term2
# })
# ll <- rowSums(ll.list)
# D <- -2 * (ll - ll.sat)


D_dR_sets <- lapply(1:length(d), function(j) {
  Dj <- d_idx2 == j
  Dj <- Dj[tOrig] & y[,2]
  dRj <- y[,1] == y[Dj, 1][1]
  list("Dj" = as.integer(which(Dj)),
       "dRj" = as.integer(which(dRj)),
       "dj" = d[j])
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

# nv <- 5
# nl <- 3
# X <- matrix(rnorm(7 * nv), nrow = 7)
# beta <- matrix(rnorm(nv * nl), nrow = nv)
# 
# idx1 <- c(1, 3)
# idx2 <- c(2, 4, 5)
# idx <- c(idx1, idx2)
# colSums(exp(X[idx1,] %*% beta)) + colSums(exp(X[idx2,] %*% beta))
# colSums(exp(X[idx,] %*% beta))

D <- -2 * (ll - ll.sat) 
D.gn <- glmnet::coxnet.deviance(y = y, x = X, beta = fit.bl$beta)
D.gn/D



















