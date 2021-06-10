library(glmnet)
library(biglasso)
library(survival)

#=======================================#
#================ SETUP ================#
#=======================================#
### GENERATE DATA 
set.seed(124)
nobs  <- 1000
nvars <- 50
p_nz  <- 0.25
beta_true <- rnorm(nvars, 0, 0.1) * rbinom(nvars, 1, p_nz)/10

dat <- coxed::sim.survdata(N = nobs, xvars = nvars, beta = beta_true)
## note that this is non-sparse in it's current form, but this should be fine for our present goals
#dat <- coxed::sim.survdata(N = nobs, xvars = nvars)

# extract data & copy it into *full incase I decide to manipulate X and y later
X <- Xfull <- as.matrix(dat$xdata)
y <- yfull <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- Xbigfull <- as.big.matrix(X)

### MODEL PARAMETERS
alpha <- 0.5
penalty <- "enet" 
lambda  <- exp(seq(-1, -6, length.out = 8)) #c(0.1, 0.01)

##### FIT MODELS
fit.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
fit.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)

beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- as.matrix(fit.gn$beta)

#beta <- beta.gn # used later when I try to reconstruct the deviance/likelihood values  

#=============================================#
#================ likelihoods ================#
#=============================================#
calc.ll0 <- function(X, y, beta) {
  # NAIVE COX LIKELIHOOD CALCULATIONS (assumes no tied times => unique sorting exists)
  #    L(beta) = prod_{i: C_i = 1} exp(X_i * beta) / (sum_{j: Y_j >= Y_i} exp(X_j * beta))
  # where C_i = 1 if the data is uncensored and C_i = 0 if the data is right-censored
  # and {j: Y_j >= Y_i} is the set of observations for which the failure/survival time
  # occurred at a time later than time Y_i. Hence,
  #  ll(beta) = sum_{i: C_i = 1} [ X_i * beta - log(sum_{j: Y_j >= Y_i} exp(X_j * beta)) ]
  
  ### saturated log-likelihood
  # code essentially copied from glmnet's implementation (see glmnet::coxnet.deviance)
  # I think this section can be written more succinctly by
  # d <- table(y[y[,2] == 1, 1])
  # ll.sat <- -sum(d * log(d))
  weights <- rep(1, nrow(y)) # currently unimplemented here, but glmnet permits obs.-wise weights
  wd <- weights[y[,2] == 1]  # weights associated with uncensored obs.
  tyd <- y[y[,2] == 1, 1]    # times associated with uncensored obs.
  if (any(duplicated(tyd))) {
    wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
  }
  wd <- wd[wd > 0]
  ll.sat <- -sum(wd * log(wd))
  
  ### fitted log-likelihood
  ll <- 0
  for (i in 1:nrow(y)) {
    if (y[i,2] == 1) {
      j <- y[,1] >= y[i,1]
      
      ll <- ll + X[i,] %*% beta - log(colSums(exp(X[j,] %*% beta)))
    }
  }
  
  ### null log-likelihood (log-likelihood corresp. to the 0 model, i.e. all coefs = 0)
  ll.null <- 0
  for (i in 1:nrow(y)) {
    if (y[i,2] == 1) {
      j <- y[,1] >= y[i,1]

      # clearly this can be written more simply, but for exposition it becomes trivially
      # apparent what the null model is doing:
      #ll.null <- ll.null + X[i,] %*% rep(0, ncol(X)) - log(sum(exp(X[j,,drop=F] %*% rep(0, ncol(X)))))
      ll.null <- ll.null - log(sum(j))
    }
  }
  # can also be written in R via the one-liner
  # ll.null <- -sum(y[,2] * log(apply(y, 1, function(yi) sum(y[,1] >= yi[1]))))
  
  return (list("ll_sat" = ll.sat, "ll_null" = ll.null, "ll" = as.numeric(ll)))
}
calc.ll1 <- function(X, y, beta) {
  ###
  ### NOTE: Should R_j (see below) include censored times or exclude them?
  ###
  # Iterative implementation of the Breslow approximation to the partial log-likelihood.
  #
  # For J unique uncensored survival/failure times, consider the j-th unique time. Let
  # D_j denote the set of (uncensored) observations with associated with j-th unique time,
  # let R_j denote the set of (uncensored) observations with failure times >= the j-th
  # unique time (observations at-risk at the j-th time). Define S_j as the sum of
  # covariates over D_j (element-wise sum),
  #   S_j = sum_{i in D_j} X_i,
  # and let d_j be the number of (uncensored) observations with the j-th unique 
  # failure time. Then,
  #   ll(beta) = sum_{j = 1 to J} [ S_j * beta - d_j * log(sum_{k in R_j} exp(X_k * beta)) ]
  #
  
  d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)
  dtime <- sort(unique(y[y[,2]==1,1]))  # sorted unique failure times for all uncensored obs.
  
  # saturated
  ll.sat <- -sum(d * log(d))
  
  # fitted
  ll <- rep(0, ncol(beta))
  for (j in 1:length(d)) {
    D_j <- (y[,1] == dtime[j]) & (y[,2] == 1)
    R_j <- (y[,1] >= dtime[j]) #& (y[,2] == 1) ### NOTE: Should this include censored?
    S <- colSums(X[D_j,,drop=F])
    ll <- ll + S %*% beta - d[j] * log(colSums(exp(X[R_j,,drop=F] %*% beta)))
  }
  
  # null
  ll.null <- 0
  for (j in 1:length(d)) {
    R_j <- (y[,1] >= dtime[j]) #& (y[,2] == 1) ### NOTE: Should this include censored?
    #ll.null <- ll.null - d[j] * log(sum(exp(X[R_j,] %*% rep(0, ncol(X)))))
    # equivalent to
    ll.null <- ll.null - d[j] * log(sum(R_j))
  }
  
  return (list("ll_sat" = ll.sat, "ll_null" = ll.null, "ll" = as.numeric(ll)))
}
calc.ll2 <- function(X, y, beta) {
  # Iterative implementation of the Breslow approximation to the partial log-likelihood,
  # retooled from the code used in biglasso (although upon reflection there must surely
  # be a more elegant way to accomplish what we're doing here).
  
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
  
  # saturated log-likelihood
  ll.sat <- -sum(d * log(d))
  
  # fitted log-likelihood
  ll <- 0
  for (j in 1:length(d)) {
    Dj <- d_idx2 == j # indices of the death/failures corresponding to the j-th unique failure time
    Dj <- Dj[tOrig] & y[,2] # hacky solution to remove censored obs & place the indices in the original ordering
    term1 <- colSums(X[Dj,,drop=F]) %*% beta
    Rj <- (y[,1] >= y[Dj,1][1]) # all times in Dj will be the same, so we pick the first
    #term2 <- d[j] * log(sum(exp(X[Rj,,drop=F] %*% beta)))
    XRbeta <- X[Rj,,drop=F] %*% beta
    term2 <- d[j] * log(colSums(exp(XRbeta)))
    
    #print(rbind(term1, term2))
    
    ll <- ll + (term1 - term2)
  }
  
  # null log-likelihood (log-likelihood corresp. to the 0 model, i.e. all coefs = 0)
  ll.null <- 0
  for (j in 1:nrow(y)) {
    if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2] (no ties)
      Rj <- (y[,1] >= y[j,1])
      term2 <- log(sum(Rj))
      ll.null <- ll.null - term2
    }
  }
  
  return (list("ll_sat" = ll.sat, "ll_null" = ll.null, "ll" = as.numeric(ll)))
}
calc.ll3 <- function(X, y, beta) {
  # REVERSE ITERATIVE CALCULATING OVER THE SET-DIFFERENCE OF THE AT-RISK SET
  # (cumulative calculation for the second term in the fitted ll)
  
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

  # saturated  
  ll.sat <- -sum(d * log(d))

  # fitted
  D_dR_sets <- lapply(1:length(d), function(j) {
    Dj <- d_idx2 == j
    Dj <- Dj[tOrig] & y[,2]
    dRj <- y[,1] == y[Dj, 1][1]
    list("Dj" = as.integer(which(Dj)),
         "dRj" = as.integer(which(dRj)),
         "dj" = d[j])
  })
  
  expXRbeta_colsum <- rep(1, ncol(beta)) # rep(0, length(lambda)) # also, rep(0, ncol(beta)) works
  ll <- 0
  for (j in length(d):1) {
    Dj <- D_dR_sets[[j]]$Dj
    dRj <- D_dR_sets[[j]]$dRj
    term1 <- colSums(X[Dj,,drop=F]) %*% beta
    expXRbeta_colsum <- expXRbeta_colsum + colSums(exp(X[dRj,,drop=F] %*% beta))
    term2 <- d[j] * log(expXRbeta_colsum)
    
    print(rbind(term1, term2))
    
    ll <- ll + (term1 - term2)
  }
  
  # null
  expXRbeta_colsum <- 0 # rep(0, length(lambda)) # also, rep(0, ncol(beta)) works
  ll.null <- 0
  for (j in length(d):1) {
    Dj <- D_dR_sets[[j]]$Dj
    dRj <- D_dR_sets[[j]]$dRj
    ###############################################
    #dRj <- setdiff(dRj, which(y[,2] == 0)) ##### TEMP
    ###############################################
    term1 <- colSums(X[Dj,,drop=F]) %*% rep(0, ncol(X))
    expXRbeta_colsum <- expXRbeta_colsum + colSums(exp(X[dRj,,drop=F] %*% rep(0, ncol(X))))
    term2 <- d[j] * log(expXRbeta_colsum)
    ll.null <- ll.null + (term1 - term2)
  }
  ll.null <- as.numeric(ll.null)
  
  return (list("ll_sat" = ll.sat, "ll_null" = ll.null, "ll" = as.numeric(ll)))
}
calc.ll.gn <- function(X, y, beta) {
  # extracting likelihoods from glmnet's deviance calculations
  dev.gn     <- glmnet::coxnet.deviance(x = X, y = y, beta = beta) # 2 * (ll_sat - ll)
  nulldev.gn <- glmnet::coxnet.deviance(x = X, y = y, beta = 0)       # 2 * (ll_sat - ll_null)
  
  ### saturated log-likelihood
  # code essentially copied from glmnet's implementation (see glmnet::coxnet.deviance)
  weights <- rep(1, nrow(y)) # currently unimplemented here, but glmnet permits obs.-wise weights
  wd <- weights[y[,2] == 1]  # weights associated with uncensored obs.
  tyd <- y[y[,2] == 1, 1]    # times associated with uncensored obs.
  if (any(duplicated(tyd))) {
    wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
  }
  wd <- wd[wd > 0]
  ll.sat <- -sum(wd * log(wd))
  
  ### fitted
  ll <- ll.sat - dev.gn/2
  
  ### null 
  ll.null <- ll.sat - nulldev.gn/2
  
  return (list("ll_sat" = ll.sat, "ll_null" = ll.null, "ll" = as.numeric(ll)))
}

# one cens, two uncens (diff times, one same as cens)
idx10 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] != 1 & y[,2] == 1)[1])
# two cens (diff times, one same as the uncens), one uncens
idx11 <- c(which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1])

idx <- idx11

ll0 <- calc.ll0(X[idx,], y[idx,], beta.bl)
ll1 <- calc.ll1(X[idx,], y[idx,], beta.bl)
ll2 <- calc.ll2(X[idx,], y[idx,], beta.bl)
ll3 <- calc.ll3(X[idx,], y[idx,], beta.bl)
ll.gn <- calc.ll.gn(X[idx,], y[idx,], beta.bl)
# ll0$ll
# ll1$ll
# ll.gn$ll

D <- rbind(-2*(ll0$ll - ll0$ll_sat),
           -2*(ll1$ll - ll1$ll_sat),
           -2*(ll2$ll - ll2$ll_sat),
           -2*(ll3$ll - ll3$ll_sat),
           -2*(ll.gn$ll - ll.gn$ll_sat),
           cox.deviance(X = Xbig, y = y, beta = fit.bl$beta, row.idx = idx)$dev)
rownames(D) <- c("ll0", "ll1", "ll2", "ll3", "ll.gn", "cox.dev")

# ll0$ll_null
# ll1$ll_null
# ll2$ll_null
# ll3$ll_null
# ll.gn$ll_null
D

ll0$ll
ll1$ll
ll2$ll
ll3$ll
ll.gn$ll


