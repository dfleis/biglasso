#####
# SIMULATED DATA
#
# Try and come up with some way to validate the cv.biglasso output
# (compute & compare with estimates provided by cv.glmnet?)
#
#
####
#=======================================#
#====== DATA GENERATING MECHANISM ======#
#=======================================#
sim.surv.weib <- function(n, lambda=0.01, rho=1, beta, rate.cens=0.001) {
  # generate survival data (Weibull baseline hazard), adapted from
  # https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
  p <- length(beta)
  X <- matrix(rnorm(n * p), nrow = n)
  
  # Weibull latent event times
  v <- runif(n = n)
  latent.time <- (-log(v)/(lambda * exp(X %*% beta)))^(1/rho)
  
  # censoring times
  cens.time <- rexp(n = n, rate = rate.cens)
  
  # follow-up times and event indicators
  time <- pmin(latent.time, cens.time)
  status <- as.numeric(latent.time < cens.time)
  
  y <- cbind(time, status)
  colnames(y) <- c("time", "status")
  
  # data set
  return (list(X = X, y = y))
}
sim.surv <- function(n, p, p_nz, seed, ...) {
  # n = number of obs.
  # p = number of covariates.
  # p_nz = proportion of non-zero coefficients (induce sparsity)
  if (!missing(seed)) set.seed(seed)
  beta <- rnorm(p, 0, 1) * rbinom(p, 1, p_nz)
  dat <- sim.surv.weib(n = n, beta = beta, ...)
  list("beta" = beta, "y" = dat$y, "X" = dat$X)
}
#=======================================#
#================ SETUP ================#
#=======================================#
library(glmnet)
library(biglasso)

### model parameters
penalty <- "enet"
alpha   <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
lambda  <- exp(seq(0, -6, length.out = 100))
ncores  <- 6

### generate data
set.seed(124)
n <- 100
p <- 1000
p_nz <- 0.1

beta <- rnorm(p, 0, 1) * rbinom(p, 1, p_nz)

nsims <- 10

#if (parallel) {doMC::registerDoMC(cores = ncores)}
pt <- proc.time()
sim <- replicate(nsims, {
  dat <- sim.surv.weib(n = n, beta = beta)
  y <- dat$y
  X <- dat$X
  Xbig <- as.big.matrix(dat$X)

  fit.bl <- biglasso(X       = Xbig,
                     y       = y,
                     family  = "cox",
                     penalty = penalty,
                     alpha   = alpha,
                     lambda  = lambda,
                     ncores  = ncores)

  fit.gn <- glmnet(x        = X,
                   y        = y,
                   family   = "cox",
                   alpha    = alpha,
                   lambda   = lambda)
                   #parallel = parallel) # I think this only exists within cv.glmnet
  
  beta.err <- sum((fit.bl$beta - fit.gn$beta)^2) # Frob norm

  beta.l2.bl <- apply(fit.bl$beta, 2, function(b) sqrt(sum((b - beta)^2)))
  beta.l1.bl <- apply(fit.bl$beta, 2, function(b) sum(abs(b - beta)))
  beta.l2.gn <- apply(fit.gn$beta, 2, function(b) sqrt(sum((b - beta)^2)))
  beta.l1.gn <- apply(fit.gn$beta, 2, function(b) sum(abs(b - beta)))
  
  nz.bl <- (fit.bl$beta != 0)
  nz.gn <- (fit.gn$beta != 0)
  
  #nz.beta.rat <- as.vector(fit.bl$beta[nz.bl & nz.gn])/as.vector(fit.gn$beta[nz.bl & nz.gn])
  
  nz.tot.bl <- apply(nz.bl, 2, sum)
  nz.tot.gn <- apply(nz.gn, 2, sum)
  
  nz.diff <- xor(nz.bl, nz.gn) # check if covariate is selected by one method but not the other
  nz.diff.lam <- apply(nz.diff, 2, sum)
  nz.diff.var <- apply(nz.diff, 1, sum)
  
  return (list("beta.err"    = beta.err, 
               "nz.tot.bl"   = nz.tot.bl, 
               "nz.tot.gn"   = nz.tot.gn, 
               "nz.diff.lam" = nz.diff.lam, 
               "nz.diff.var" = nz.diff.var,
               "beta.l2.bl"  = beta.l2.bl,
               "beta.l1.bl"  = beta.l1.bl,
               "beta.l2.gn"  = beta.l2.gn,
               "beta.l1.gn"  = beta.l1.gn))
})
(tm <- proc.time() - pt)

beta.err <- unlist(sim["beta.err",])
nz.tot.bl <- simplify2array(sim["nz.tot.bl",])
nz.tot.gn <- simplify2array(sim["nz.tot.gn",])
nz.diff.lam <- simplify2array(sim["nz.diff.lam",])
nz.diff.var <- simplify2array(sim["nz.diff.var",])
beta.l2.bl <- simplify2array(sim["beta.l2.bl",])
beta.l1.bl <- simplify2array(sim["beta.l1.bl",])
beta.l2.gn <- simplify2array(sim["beta.l2.gn",])
beta.l1.gn <- simplify2array(sim["beta.l1.gn",])


myclrs <- c(rgb(26/255, 133/255, 255/255, 0.7), 
            rgb(212/255, 17/255, 89/255, 0.7))
plot(NA, xlim = range(lambda), ylim = range(nz.tot.bl, nz.tot.gn), log = 'x',
     xlab = expression(lambda), ylab = "Nb. Vars.", 
     main = "Number of Variables Retained by\nPenalized Cox Implementations")
grid(); abline(h = 0, lwd = 1.5, col = 'gray50')
for (i in 1:nsims) {
  lines(nz.tot.bl[,i] ~ lambda, type = 's', col = myclrs[1], lwd = 2)
  lines(nz.tot.gn[,i] ~ lambda, type = 's', col = myclrs[2], lwd = 2)
}
legend("topright", legend = c("biglasso", "glmnet"), seg.len = 2, col = myclrs, lty = "solid", lwd = 4)

nz.tot.diff <- nz.tot.gn - nz.tot.bl
plot(NA, xlim = range(lambda), ylim = range(nz.tot.diff), log = 'x',
     xlab = expression(lambda), ylab = "nvars(biglasso) - nvars(glmnet)", 
     main = "Difference between the Nb. of Covariates\nRetained by biglasso and glmnet")
grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray50')
for (i in 1:nsims) {
  lines(nz.tot.diff[,i] ~ lambda, type = 's', col = rgb(0, 0, 0, 0.75), lwd = 1)
}
lines(apply(nz.tot.diff, 1, mean) ~ lambda, col = 'red', lwd = 2)

plot(NA, xlim = range(lambda), ylim = range(nz.diff.lam), log = 'x',
     xlab = expression(lambda), ylab = "Nb. Vars.", 
     main = "Number of Variables Retained by\nOne Method but not the Other")
grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray50')
for (i in 1:nsims) {
  lines(nz.diff.lam[,i] ~ lambda, type = 's', col = rgb(0, 0, 0, 0.75), lwd = 1)
}
lines(apply(nz.diff.lam, 1, mean) ~ lambda, col = 'red', lwd = 2)

# plot(NA, xlim = range(lambda), ylim = range(beta.l2.bl, beta.l2.gn), log = 'x',
#      xlab = expression(lambda), ylab = "",
#      main = "Fitted Coefficient Loss\nl_2-Norm")
# title(ylab = expression("\u2016"*beta - hat(beta)*"\u2016"[2]), line = 2.5)
# grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray50')
# for (i in 1:nsims) {
#   lines(beta.l2.bl[,i] ~ lambda, type = 'l', col = myclrs[1], lwd = 1)
#   lines(beta.l2.gn[,i] ~ lambda, type = 'l', col = myclrs[2], lwd = 1)
# }
# legend("topright", legend = c("biglasso", "glmnet"), seg.len = 2, col = myclrs, lty = "solid", lwd = 4)
# 
# plot(NA, xlim = range(lambda), ylim = range(beta.l1.bl, beta.l1.gn), log = 'x',
#      xlab = expression(lambda), ylab = "",
#      main = "Fitted Coefficient Loss\nl_1-Norm")
# title(ylab = expression("\u2016"*beta - hat(beta)*"\u2016"[1]), line = 2.5)
# grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray50')
# for (i in 1:nsims) {
#   lines(beta.l1.bl[,i] ~ lambda, type = 'l', col = myclrs[1], lwd = 1)
#   lines(beta.l1.gn[,i] ~ lambda, type = 'l', col = myclrs[2], lwd = 1)
# }
# legend("topright", legend = c("biglasso", "glmnet"), seg.len = 2, col = myclrs, lty = "solid", lwd = 4)



















