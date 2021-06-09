library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)

################### GENERATE DATA 
set.seed(124)
no <- 2500
nx <- 1000
dat <- coxed::sim.survdata(N = no, xvars = nx)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

# model & cv parameters
penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(0, -6, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- T # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 5
trace.it <- 1

################### RUN CV MODELS

# lapply(list.files("../cv-biglasso-cox/R/", full.names = T), source)
# 
# pt <- proc.time()
# NOTE: This one has to be run first if we want to keep the cv.ind/foldid the same 
# across the other cv functions (I didn't code in a way to input cv fold indices
# in this earlier version)
# cv.blR <- cv.biglasso.cox(x        = Xbig,
#                           y        = y,
#                           penalty  = penalty,
#                           alpha    = alpha,
#                           lambda   = lambda,
#                           nfolds   = nfolds,
#                           grouped  = grouped,
#                           parallel = parallel)
# proc.time() - pt

pt <- proc.time()
cv.bl <- cv.biglasso(X       = Xbig,                               
                     y       = y,
                     family  = "cox",
                     penalty = penalty,
                     alpha   = alpha,
                     nfolds  = nfolds,
                     lambda  = lambda,
                     grouped = grouped,
                     ncores  = ncores,
                     trace   = as.logical(trace.it))
proc.time() - pt


if (parallel) {doMC::registerDoMC(cores = ncores)}

pt <- proc.time()
cv.gn <- cv.glmnet(x        = X,                               
                   y        = y,
                   family   = "cox",
                   alpha    = alpha,
                   nfolds   = nfolds,
                   lambda   = lambda,
                   grouped  = grouped,
                   parallel = parallel,
                   trace.it = trace.it,
                   foldid   = cv.bl$cv.ind) 
proc.time() - pt




#############################

plot(cv.bl)
plot(cv.gn)


plot(cv.bl$cve ~ lambda, log = 'x', type = 'l', lwd = 2)
lines(cv.gn$cvm ~ lambda, col = 'red', lty = 'dashed', lwd = 2)

#plot(cv.gn$cvm - cv.bl$cve ~ lambda, log = 'x', type = 'l', lwd = 2)
