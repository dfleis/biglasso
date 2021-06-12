###############################################################
# SIMULATED DATA 
#
# Comparing glmnet & the (proposed) biglasso CV outputs
# for penalized Cox models
#
###############################################################
#=======================================#
#====== DATA GENERATING MECHANISM ======#
#=======================================#
sim.surv.weib <- function(n, lambda=0.01, rho=1, beta, rate.cens=0.001, round.times=F, round.digits=0) {
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
  if (round.times) time <- floor(time * 10^round.digits)/(10^round.digits) + 1#round(time, round.digits)
  status <- as.numeric(latent.time < cens.time)
  
  y <- cbind(time, status)
  colnames(y) <- c("time", "status")
  
  # data set
  return (list(X = X, y = y))
}
sim.surv <- function(n, p, p_nz, seed) {
  if (!missing(seed)) set.seed(seed)
  
  beta <- rnorm(p, 0, 1) * rbinom(p, 1, p_nz)
  dat <- sim.surv.weib(n = n, beta = beta)
  list("beta" = beta, "y" = dat$y, "X" = dat$X)
}
#=======================================#
#================ SETUP ================#
#=======================================#
#remove.packages("biglasso")
#devtools::install_github("dfleis/biglasso")
library(biglasso)

set.seed(124)
n <- 100
p <- 500
p_nz <- 0.1
beta <- rnorm(p, 0, 1) * rbinom(p, 1, p_nz)

dat <- sim.surv.weib(n = n, beta = beta, rho = 10, rate.cens = 0.05, round.times=T)
y <- dat$y
X <- dat$X
Xbig <- as.big.matrix(X)

table(beta != 0)
table(y[,1], y[,2])
table(y[,2])

#===========================================#
#================ RUN TESTS ================#
#===========================================#
set.seed(124) # necessary for testing as the bigmemory package has an (unintended?) effect on random number generation
penalty  <- "enet"
alpha    <- 0.1 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 10
lambda   <- exp(seq(1, -5, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
ncores   <- 1
trace.it <- 1
foldid   <- sample(cut(1:nrow(y), breaks = nfolds, labels = F))

# fit <- biglasso::biglasso(X = Xbig, y = y, family = "cox", penalty = "enet", alpha = alpha, lambda = lambda)
# biglasso::cox.deviance(X = Xbig, y = y, beta = fit$beta, row.idx = 1:nrow(y))


pt <- proc.time()
cv.bl <- biglasso::cv.biglasso(
                     X       = Xbig,                               
                     y       = y,
                     family  = "cox",
                     penalty = penalty,
                     alpha   = alpha,
                     nfolds  = nfolds,
                     lambda  = lambda,
                     grouped = grouped,
                     cv.ind  = foldid,
                     ncores  = ncores,
                     trace   = as.logical(trace.it))
proc.time() - pt



plot(cv.bl)
