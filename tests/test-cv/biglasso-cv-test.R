library(bigmemory)
library(glmnet)
library(biglasso)
library(survival)

set.seed(124)
### GENERATE DATA 
no <- 100
nx <- 25
dat <- coxed::sim.survdata(N = no, xvars = nx)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

# model & cv parameters
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(-1, -7, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- F # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
trace.it <- 1

### RUN CROSS-VALIDATED MODELS
pt <- proc.time()
cv.bl <- cv.biglasso(X      = Xbig,                               
                     y      = y,
                     family = "cox",
                     alpha  = alpha,
                     nfolds = nfolds,
                     lambda = lambda,
                     ncores = 1,
                     trace  = as.logical(trace.it))
proc.time() - pt

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


plot(cv.bl)
plot(cv.gn)
