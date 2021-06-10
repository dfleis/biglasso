##
# I noticed a bug (?) (June 9, 2021) in the deviance calculation in which
# its outputs will appreciably differ from those of glmnet::coxnet.deviance
# when a unique survival time occurs right-censored (y[obs. nb., 2] == 0).
# That is, for some right-censored survival time t, there is no uncensored
# observation with the same survival time.
#
# It looks like a single uncensored obs. (with the same time) is 
# sufficient to make the cox.deviance output identical with that of 
# glmnet::coxnet.deviance
#
# NOTES from glmnet:
# "if death times are tied with censored times, we assume the censored times occurred 
#  just before the death times in computing the Breslow approximation; if users prefer 
#  the usual convention of after, they can add a small number to all censoring times to 
#  achieve this effect."
#
#
###
#=======================================#
#================ SETUP ================#
#=======================================#
#remove.packages("biglasso")
#devtools::install_github("dfleis/biglasso")
library(biglasso)
library(glmnet)

set.seed(124)
n <- 1000
p <- 100
p_nz <- 0.5
beta <- rnorm(p, 0, 0.1) * rbinom(p, 1, p_nz)/10

dat <- coxed::sim.survdata(N = n, T = 100, xvars = p, beta = beta)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

#===========================================#
#================ RUN TESTS ================#
#===========================================#
### I don't think I need to run full CV's for the deviance tests below
### Instead, I think just getting coefficients from biglasso and glmnet
### would be fine
penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(-2, -6, length.out = 100))
grouped  <- F # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- F # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 1
trace.it <- 1

set.seed(124) # necessary for testing as the bigmemory package has an (unintended?) effect on random number generation
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

############## compare deviance calculations
set.seed(124)
idx1 <- sort(sample(1:nrow(y), 17))
set.seed(124)
idx2 <- sort(sample(1:nrow(y), 18))

idx <- idx1

dev.bl0_raw <- biglasso::cox.deviance(X = Xbig, y = y, beta = cv.bl$fit$beta, row.idx = idx)
dev.bl1_raw <- biglasso::cox.deviance(X = Xbig, y = y, beta = cv.gn$glmnet.fit$beta, row.idx = idx)
dev.bl0 <- dev.bl0_raw$dev
dev.bl1 <- dev.bl1_raw$dev
dev.gn0 <- glmnet::coxnet.deviance(y = y[idx,], x = X[idx,], beta = cv.bl$fit$beta)
dev.gn1 <- glmnet::coxnet.deviance(y = y[idx,], x = X[idx,], beta = cv.gn$glmnet.fit$beta)

plot(NA, xlab = "lambda", ylab = "Dev", log = 'x',
     xlim = range(lambda), ylim = range(dev.bl0, dev.bl1, dev.gn0, dev.gn1))
grid(nx = 5, ny = 5)
lines(dev.bl0 ~ lambda, lwd = 2)
lines(dev.bl1 ~ lambda, lwd = 2, col = 'magenta')
lines(dev.gn0 ~ lambda, lwd = 2, lty = 'dashed')
lines(dev.gn1 ~ lambda, lwd = 2, col = 'magenta', lty = 'dashed')


############## compare deviance calculations
# set.seed(124)
# idx1 <- sort(sample(1:nrow(y), 17))
# set.seed(124)
# idx2 <- sort(sample(1:nrow(y), 18))
# 
# idx3 <- c(which(y[,2] == 0)[1], which(y[,2] == 1)[1])[1]
# idx4 <- c(which(y[,1] == 96 & y[,2] == 0), which(y[,1] != 96 & y[,2] == 1)[1:2])
#
# table(y[idx1,1], y[idx1,2])
# table(y[idx2,1], y[idx2,2])
# table(y[idx3,1], y[idx3,2])
# table(y[idx4,1], y[idx4,2])
# table(y[idx5,1], y[idx5,2])

table(y[,1],y[,2])
idx1 <- which(y[,1] == 1 & y[,2] == 1)[1:2]                                     # two uncensored obs. same times
idx2 <- c(which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] != 1 & y[,2] == 1)[1]) # two uncensored obs. diff times
idx3 <- which(y[,1] == 1 & y[,2] == 0)[1:2]                                     # two censored obs. same times (should fail in both glmnet and biglasso)
idx4 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1]) # two censored obs. diff times (should fail in both glmnet and biglasso)
idx5 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] == 1 & y[,2] == 1)[1]) # one cens. one uncens, same times
idx6 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 1)[1]) # one cens. one uncens, diff times

idx7 <- c(which(y[,1] == 96 & y[,2] == 0)[1:2], which(y[,1] == 96 & y[,2] == 1)[1]) # two censored, one uncensored, same times
idx8 <- c(which(y[,1] == 1 & y[,2] == 0)[1:2], which(y[,1] != 1 & y[,2] == 1)[1]) # one cens. one uncens, diff times

# two cens (diff times), one uncens (same time as one cens)
idx9 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 1)[1])

# one cens, two uncens (one same time as cens, one diff) 
idx10 <- c(which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] != 1 & y[,2] == 1)[1])

# two cens (diff times, one same as the uncens), one uncens
idx11 <- c(which(y[,1] == 1 & y[,2] == 1)[1], which(y[,1] == 1 & y[,2] == 0)[1], which(y[,1] != 1 & y[,2] == 0)[1])


## notes: idx10 should lead to the same result for biglasso and glmnet deviance
##        idx11 should lead to different results (a time associated with a censored obs. is not shared with an uncens. obs)

y[idx1,]
y[idx2,]
y[idx3,]
y[idx4,]
y[idx5,]
y[idx6,]
y[idx7,]
y[idx8,]
y[idx9,]
y[idx10,]
y[idx11,]
idx <- idx11

dev.bl0_raw <- biglasso::cox.deviance(X = Xbig, y = y, beta = cv.bl$fit$beta, row.idx = idx)
dev.bl1_raw <- biglasso::cox.deviance(X = Xbig, y = y, beta = cv.gn$glmnet.fit$beta, row.idx = idx)
dev.bl0 <- dev.bl0_raw$dev
dev.bl1 <- dev.bl1_raw$dev
dev.gn0 <- glmnet::coxnet.deviance(y = y[idx,], x = X[idx,], beta = cv.bl$fit$beta)
dev.gn1 <- glmnet::coxnet.deviance(y = y[idx,], x = X[idx,], beta = cv.gn$glmnet.fit$beta)


# #-2 * (dev.bl0_raw$loglik - dev.bl0_raw$loglik_sat)
#dev.bl0 - dev.gn0
dev.bl0
dev.gn0
# 
# -2 * dev.bl0_raw$loglik
# dev.bl0

plot(NA, xlab = "lambda", ylab = "Dev", log = 'x',
     xlim = range(lambda), ylim = range(dev.bl0, dev.bl1, dev.gn0, dev.gn1))
grid(nx = 5, ny = 5)
lines(dev.bl0 ~ lambda, lwd = 2)
lines(dev.bl1 ~ lambda, lwd = 2, col = 'magenta')
lines(dev.gn0 ~ lambda, lwd = 2, lty = 'dashed')
lines(dev.gn1 ~ lambda, lwd = 2, col = 'magenta', lty = 'dashed')

############
# ty <- y[,"time"]
# tevent <- y[,"status"]
# ty <- ty + (1 - tevent) * 100 * .Machine$double.eps #









