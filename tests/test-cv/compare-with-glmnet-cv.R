#####
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

### cv parameters
penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(0, -6, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
ncores   <- 1

### generate data
seed <- 124 
set.seed(seed) # we must call set.seed more than once since using as.big.matrix (unintentionally?) affects RNG
Ns <- floor(10^seq(1.5, 3, length.out = 20))
Ps <- floor(10^seq(1.5, 3, length.out = 20))
p_nz <- 1

parallel <- TRUE
ncores <- 6
if (parallel) {doMC::registerDoMC(cores = ncores)}

CV.ERR <- LAM.ERR <- matrix(NA, nrow = length(Ns), ncol = length(Ps))

pto <- proc.time()
for (i in 1:length(Ns)) {
  n <- Ns[i]
  foldid <- sample(cut(1:n, breaks = nfolds, labels = F))
  
  for (j in 1:length(Ps)) {
    p <- Ps[j]
    print(sprintf("(n, p) = (%i, %i)", n, p))
    print(sprintf("(i, k) = (%i, %i)", i, j))
    
    
    print("Generating data...")
    dat <- sim.surv(n, p, p_nz)
    y <- dat$y
    X <- dat$X
    Xbig <- as.big.matrix(dat$X)
    
    print("Running cv.biglasso...")
    pt <- proc.time()
    cv.bl <- cv.biglasso(X       = Xbig,
                         y       = y,
                         family  = "cox",
                         penalty = penalty,
                         alpha   = alpha,
                         nfolds  = nfolds,
                         lambda  = lambda,
                         grouped = grouped,
                         cv.ind  = foldid,
                         ncores  = ncores)
    print(proc.time() - pt)
    
    print("Running cv.glmnet...")
    pt <- proc.time()
    cv.gn <- cv.glmnet(x        = X,
                       y        = y,
                       family   = "cox",
                       alpha    = alpha,
                       nfolds   = nfolds,
                       lambda   = lambda,
                       grouped  = grouped,
                       parallel = parallel,
                       foldid   = foldid)
    print(proc.time() - pt)
    
    cv.err.pct <- 1 - cv.gn$cvm/cv.bl$cve
    max.idx <- which.max(abs(cv.err.pct))
    LAM.ERR[i,j] <- lambda[max.idx]
    CV.ERR[i,j]  <- cv.err.pct[max.idx]
  }
}
proc.time() - pto

#fields::image.plot(Ns, Ps, log(LAM.ERR), col = pals::parula(33), zlim = range(log(lambda)))
#fields::image.plot(Ns, Ps, CV.ERR, col = pals::parula(33), log = 'xy')
fields::image.plot(log10(Ns), log10(Ps), log(LAM.ERR), col = pals::parula(33), zlim = range(log(lambda)))
fields::image.plot(log10(Ns), log10(Ps), CV.ERR, 
                   col = pals::coolwarm(33), 
                   zlim = c(-max(abs(CV.ERR)), max(abs(CV.ERR))))









