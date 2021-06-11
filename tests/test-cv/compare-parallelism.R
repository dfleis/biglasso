#####
# Verify that the biglasso outputs are identical when using
# multithreading (via 'ncores' argument)
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
  time <- round(pmin(latent.time, cens.time)) + 1
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

set.seed(124)
### data parameters
n <- 1000
p <- 5
p_nz <- 1

### biglasso parameters
penalty <- "enet"
alpha   <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
lambda  <- exp(seq(0, -6, length.out = 100))
nfolds  <- 5
grouped <- T
foldid  <- sample(cut(1:n, breaks = nfolds, labels = F))
ncores  <- 4

### generate data
dat <- sim.surv(n = n, p = n, p_nz = p_nz, rho = 10, rate.cens=0.1)
y <- dat$y
X <- dat$X
Xbig <- as.big.matrix(X)
table(y[,1], y[,2])
table(y[,2])

pt <- proc.time()
cv.bl0 <- cv.biglasso(X       = Xbig,
                      y       = y,
                      family  = "cox",
                      penalty = penalty,
                      alpha   = alpha,
                      lambda  = lambda,
                      nfolds  = nfolds,
                      grouped = grouped,
                      cv.ind  = foldid,
                      ncores  = 1,
                      trace   = T)
proc.time() - pt
pt <- proc.time()
cv.bl1 <- cv.biglasso(X       = Xbig,
                      y       = y,
                      family  = "cox",
                      penalty = penalty,
                      alpha   = alpha,
                      lambda  = lambda,
                      nfolds  = nfolds,
                      grouped = grouped,
                      cv.ind  = foldid,
                      ncores  = ncores,
                      trace   = T)
proc.time() - pt

doMC::registerDoMC(cores = ncores)
pt <- proc.time()
cv.gn <- cv.glmnet(x        = X,
                   y        = y,
                   family   = "cox",
                   alpha    = alpha,
                   nfolds   = nfolds,
                   lambda   = lambda,
                   grouped  = grouped,
                   parallel = T,
                   foldid   = foldid,
                   trace.it = 1)
proc.time() - pt


#lapply(list.files("~/projects/cv-biglasso-cox/R/", full.names = T), source)
#plot.compare.cv2(cv.bl1, cv.gn)
plot.compare.cv2(cv.bl0, cv.gn)
