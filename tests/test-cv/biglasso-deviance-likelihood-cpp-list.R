library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)
library(Rcpp)

set.seed(124)
### GENERATE DATA 
no <- 1000
nx <- 15
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
lambda  <- exp(seq(-1, -6, length.out = 100)) #c(0.1, 0.01)

##### FIT MODELS
#fit.cx <- coxph(Surv(y) ~ X, ties = "breslow")
fit.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
fit.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)

beta.cx <- coef(fit.cx)
beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- as.matrix(fit.gn$beta)

beta <- beta.bl # used later when I try to reconstruct the deviance/likelihood values  

##### likelihoods
d <- as.numeric(table(y[y[,2]==1,1])) # counts of unique failure times for all uncensored obs. (y[,2] == 1)

#satdev <- 2 * sum(d * log(d))
ll.sat <- -sum(d * log(d))
ll.cx  <- fit.cx$loglik

dev.cx <- -2 * (ll.cx - ll.sat)
dev.bl <- fit.bl$loss
dev.gn <- (1 - fit.gn$dev.ratio) * fit.gn$nulldev
dev.gn.null <- fit.gn$nulldev # defined as 2 * (ll_sat - ll_null)

pt <- proc.time()
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

D_R_sets <- lapply(1:length(d), function(j) { 
  # NOTE: structure is slightly different in the cox.deviance function
  # so that it's more amenable to passing it to C++, i.e.
  # indexes starting at 0, unnecessary inclusion of d[j].
  Dj <- d_idx2 == j
  Dj <- Dj[tOrig] & y[,2]
  Rj <- (y[,1] >= y[Dj,1][1])
  list("Dj" = as.integer(which(Dj)), 
       "Rj" = as.integer(which(Rj)),
       "dj" = d[j]) 
})

ll.list <- sapply(D_R_sets, function(DRs) {
  term1 <- colSums(X[DRs$Dj,,drop=F]) %*% beta
  term2 <- DRs$d * log(colSums(exp(X[DRs$Rj,,drop=F] %*% beta)))
  term1 - term2
})
ll <- rowSums(ll.list)
D <- -2 * (ll - ll.sat)
proc.time() - pt

###### SEPARATE LOGLIKELIHOOD/DEVIANCES FUNCTIONS
idx <- 1:nrow(y)

pt <- proc.time()
out <- cox.deviance(X = Xbig, y = y, beta = fit.bl$beta, row.idx = idx)
proc.time() - pt

pt <- proc.time()
out.gn <- glmnet::coxnet.deviance(y = y[idx,], x = X[idx,], beta = fit.bl$beta)
proc.time() - pt


plot(out.gn ~ lambda, log = 'x', type = 'l', lwd = 2)
lines(out$D ~ lambda, lwd = 2, col = 'red')
lines(D ~ lambda, lwd = 2, col = 'blue')
plot((out$D - out.gn)/out$D ~ lambda, log = 'x', type = 'l')
plot((out$D - D)/out$D ~ lambda, log = 'x', type = 'l')







cppFunction(depends = "RcppArmadillo", "
SEXP fx(SEXP M) {
  arma::mat m = Rcpp::as<arma::mat>(M);
  arma::mat colsum = arma::sum(m, 0);
  return List::create(colsum);
}
")
A <- matrix(rnorm(4 * 3), nrow = 4)
colSums(A)
fx(A)

Amat <- as.matrix(A)
matrix(1, ncol = nrow(Amat)) %*% Amat 
colSums(Amat)

cppFunction(plugins = "type_info", "
void printClass(SEXP x) {
  std::cout << \"hello, world!\" << std::endl; 
  //std::cout << x << std::endl;  
}
")
printClass(1)

#sourceCpp("tests/test-cpp/test-cpp.cpp")
cppFunction(depends = "RcppArmadillo", "
SEXP fx() {
  arma::mat Z = arma::zeros(3, 10);
  Z(0, 5) = 1.0;
  return List::create(Z);
}
")
fx()

cppFunction(depends = "RcppArmadillo", "
SEXP fx(SEXP M, double x) {
  arma::mat m = Rcpp::as<arma::mat>(M);
  return List::create(x * m);
}
")
X <- matrix(runif(4 * 3), nrow = 4)
3.1 * X
fx(X, 3.1)

cppFunction(depends = "RcppArmadillo", "
void fx(SEXP M) {
  arma::mat m = Rcpp::as<arma::mat>(M);
  m.print();
}
")
v <- matrix(rnorm(4 * 3), nrow = 4)
v
fx(v)

################ manual tests 
ll.sat <- -sum(d * log(d))

D_R_sets <- lapply(1:length(d), function(j) {
  Dj <- d_idx2 == j
  Dj <- Dj[tOrig] & y[,2]
  Rj <- (y[,1] >= y[Dj,1][1])
  list("Dj" = as.integer(which(Dj)), 
       "Rj" = as.integer(which(Rj)))
})


ll.list <- sapply(D_R_sets, function(DRs) {
  term1 <- colSums(X[DRs$Dj,,drop=F]) %*% beta
  term2 <- DRs$d * log(colSums(exp(X[DRs$Rj,,drop=F] %*% beta)))
  term1 - term2
})
D <- -2 * (rowSums(ll.list) - ll.sat)
  
ll.null <- 0
for (j in 1:nrow(y)) {
  if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2] (no ties)
    Rj <- (y[,1] >= y[j,1])
    term2 <- log(sum(Rj))
    ll.null <- ll.null - term2
  }
}
D0 <- -2 * (ll.null - ll.sat)
  

# pt <- proc.time()
# devs1L <- calc.devs.LIST(X, y, beta.bl)
# proc.time() - pt
# pt <- proc.time()
# devs2L <- calc.devs.LIST(X, y, beta.gn)
# proc.time() - pt
# 
# 
# devs1
# devs1L
