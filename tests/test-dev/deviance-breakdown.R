#=======================================#
#================ SETUP ================#
#=======================================#
library(biglasso)
library(glmnet)
library(survival)

############# GENERATE DATA #############
set.seed(124)
nobs  <- 1000
nvars <- 100
p_nz <- 0.5
beta <- rnorm(nvars, 0, 0.1) * rbinom(nvars, 1, p_nz)/10

dat <- coxed::sim.survdata(N = nobs, T = 100, xvars = nvars, beta = beta)
dat$data$y <- dat$data$y
X <- Xfull <- as.matrix(dat$xdata)
y <- yfull <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- Xbigfull <- as.big.matrix(X)

############# FIT MODELS #############
### model parameters
alpha <- 0.5
penalty <- "enet" 
lambda  <- c(1, 0.1, 0.01, 0.001, 0) #exp(seq(-1, -6, length.out = 7)) #c(0.1, 0.01)

### Cox regressions
fit.cx <- coxph(Surv(y) ~ X, ties = "breslow")
fit.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
fit.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)

### extract coefs
beta.cx <- coef(fit.cx)
beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- as.matrix(fit.gn$beta)

beta <- beta.bl

### select problematic subsets
# one cens, two uncens (one same time as cens, one diff) 
idx10 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] != 1 & y[,2] == 1)[1])
# two cens (diff times, one same as the uncens), one uncens
idx11 <- c(which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1])


idx <- idx11
y <- y[idx,]
X <- X[idx,]
Xbig <- as.big.matrix(X)
############# CALCULATE LIKELIHOOD & DEVIANCE #############
table(y[,"time"], y[,"status"]) # summary of failure times split by censor/uncensored status
# following two variables are used in glmnet
ty <- y[,"time"] 
tevent <- y[,"status"]

### calculate basic log-likelihoods
ll <- 0
for (i in 1:nrow(y)) {
  if (y[i,"status"] == 1) {
    j <- y[,"time"] >= y[i,"time"]
    
    ll <- ll + X[i,] %*% beta - log(colSums(exp(X[j,] %*% beta)))
  }
}

ll.null <- 0
for (i in 1:nrow(y)) {
  if (y[i,"status"] == 1) {
    j <- y[,"time"] >= y[i,"time"]

    # clearly this can be written more simply, but for exposition it's trivially obvious what the null model is doing
    #ll.null <- ll.null + X[i,] %*% rep(0, ncol(X)) - log(sum(exp(X[j,,drop=F] %*% rep(0, ncol(X))))) 
    # equivalent to
    ll.null <- ll.null - log(sum(j))
  }
}
# null likelihood can also be written as
# -sum(y[,"status"] * log(apply(y, 1, function(yi) sum(y[,"time"] >= yi["time"]))))


### extract likelihood/deviance values from the various fitted model outputs
ll.cx  <- fit.cx$loglik # coxph()
# 'loglik' = a vector of length 2 containing the log-likelihood with the initial 
# values and with the final values of the coefficients.
dev.bl <- fit.bl$loss



### GLMNET
nulldev.gn <- fit.gn$nulldev # 2 * (ll_sat - ll_null) # coxnet.deviance(y = y, x = X, beta = 0)
dev.gn <- (1 - fit.gn$dev.ratio) * nulldev.gn # 2 * (ll_sat - ll)

# ll.null.gn <- -(nulldev.gn/2 - ll.sat.gn)
# ll.gn <- -(dev.gn/2 - ll.sat.gn)
# 
# ll.sat.gn
# ll.null.gn
# ll.gn
# how glmnet calculates saturated likelihood
weights <- rep(1, nrow(y))
wd <- weights[tevent == 1]
tyd <- ty[tevent == 1]
if (any(duplicated(tyd))) {
  wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
}
wd <- wd[wd > 0]

ll.sat.gn <- -sum(wd * log(wd))

-2 * (ll - ll.sat.gn)


# idx <- idx11
cox.deviance(X = Xfull[idx,], y = yfull[idx,], beta = fit.bl$beta, row.idx = idx)$dev
coxnet.deviance(x = Xfull[idx,], y = yfull[idx,], beta = fit.bl$beta)






d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)

-sum(d * log(d))



ll.sat <- -sum(d * log(d))
ll.cx  <- fit.cx$loglik

dev.cx <- -2 * (ll.cx - ll.sat)
dev.bl <- fit.bl$loss
dev.gn <- (1 - fit.gn$dev.ratio) * fit.gn$nulldev
dev.gn.null <- fit.gn$nulldev # defined as 2 * (ll_sat - ll_null)



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
  
  return (D)
  #return (c("nulldev" = D0, "dev" = D, "ll_sat" = ll.sat, "ll" = ll, "ll_null" = ll.null))
}

# Dmat <- matrix(c("cx" = dev.cx, "bl" = c(NA, dev.bl), "gn" = c(dev.gn.null, dev.gn)), ncol = 2, nrow = 3, byrow = T)
# colnames(Dmat) <- c("D0", "D")
# rownames(Dmat) <- c("cx", "bl", "gn")

pt <- proc.time()
devs1 <- calc.devs(X, y, beta.bl)
proc.time() - pt
pt <- proc.time()
devs2 <- calc.devs(X, y, beta.gn)
proc.time() - pt

#round(devs1, 4)
#round(devs2, 4)






