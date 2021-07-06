#=======================================#
#================ SETUP ================#
#=======================================#
library(glmnet)
library(biglasso)
library(coxed)

set.seed(124)

### cv parameters
penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(2, -1, length.out = 50))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- T # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 5
trace.it <- 1

### data parameters & generate data
n <- 100
p <- 5000
p_nz <- 0.90 # proportion of nonzerom coefs.

# I run sim.survdata twice so that the coxed library generates plausible beta
# values and then I 'sparsify' the betas and re-run sim.survdata. Otherwise
# if I specify the sparse betas outright it tends to break
data <- coxed::sim.survdata(N = n, xvars = p)
nb_nz    <- ceiling(p_nz * p)
beta_ind <- rep(c(0, 1), c(p - nb_nz, nb_nz))
beta     <- data$betas * sample(beta_ind)
data <- coxed::sim.survdata(X = data$xdata, beta = beta)

y <- cbind("time" = data$data$y, "status" = data$data$failed)
X <- as.matrix(data$xdata)
Xbig <- as.big.matrix(X)


################### RUN CV MODELS
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


#######################

plot(cv.gn)
plot(cv.bl)




