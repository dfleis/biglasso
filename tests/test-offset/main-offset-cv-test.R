##############################################################################
# Testing whether the implementation of an offset argument for biglasso's
# Cox regression works as expected.
#
#
###############################################################################
library(glmnet)
library(biglasso)
library(coxed)

set.seed(124)

### data parameters & generate data
n <- 1000
p <- 150
p_nz  <- 0.1 # proportion of nonzerom coefs.
nb_nz <- ceiling(p_nz * p)

beta_ind <- sample(rep(c(0, 1), c(p - nb_nz, nb_nz)))
beta <- rnorm(p) * beta_ind

X <- matrix(rnorm(n * p), nrow = n)
eta <- X %*% beta
time <- sapply(eta, function(e) rexp(1, rate = exp(e)))
status <- rbinom(n, 1, 0.75)
y <- cbind("time" = time, "status" = status)

offset <- rnorm(n)*10


### model params
penalty <- "enet"
alpha   <- 0.5
lambda  <- exp(seq(1, -7, length.out = 100))
nfolds  <- 5
foldid  <- ceiling(sample(1:n)/n*nfolds)
grouped <- T



Xbig <- as.big.matrix(X)

# data <- coxed::sim.survdata(N = n, xvars = p, beta = beta)
# #table(data$data$y, data$data$failed)
# 
# X <- as.matrix(data$xdata)
# Xbig <- as.big.matrix(X)
# y <- cbind("time" = data$data$y, "status" = data$data$failed)
# 


#======================== FIT MODELS ========================#


### run models
cv.gn.off <- cv.glmnet(x       = X, 
                       y       = y,
                       offset  = offset,
                       family  = "cox", 
                       alpha   = alpha,
                       lambda  = lambda,
                       nfolds  = nfolds,
                       foldid  = foldid,
                       grouped = grouped)

cv.bl.off <- cv.biglasso(X       = Xbig,
                         y       = y,
                         family  = "cox",
                         offset  = offset,
                         penalty = penalty,
                         alpha   = alpha,
                         lambda  = lambda,
                         nfolds  = nfolds,
                         cv.ind  = foldid,
                         grouped = grouped)


source("~/projects/cv-biglasso-cox/R/plot.cv.biglasso.cox.R")
source("~/projects/cv-biglasso-cox/R/getmin.lambda.R")

plot.compare.cv2(cv.bl, cv.gn)
plot.compare.cv2(cv.bl.off, cv.gn.off)

beta.gn.off <- as.matrix(cv.gn.off$glmnet.fit$beta)
beta.bl.off <- as.matrix(cv.bl.off$fit$beta)

beta.diff <- beta.gn.off - beta.bl.off

err <- apply(beta.diff, 2, function(b) sqrt(sum(b^2)))
plot(err ~ log(lambda), type = 'l')


# coxnet.deviance(y = y, x = X, offset = offset, beta = beta.bl.off)
# cox.deviance(X = Xbig, y = y, offset = offset, beta = beta.bl.off)$dev

lambda <- exp(seq(0, -5, length.out = 100))

idx <- which(foldid == 1)
fit.gn <- glmnet(x       = X[idx,], 
                 y       = y[idx,],
                 offset  = offset[idx],
                 family  = "cox", 
                 alpha   = alpha,
                 lambda  = lambda)
fit.bl1 <- biglasso(X       = Xbig,
                    y       = y,
                    family  = "cox",
                    offset  = offset,
                    penalty = penalty,
                    alpha   = alpha,
                    lambda  = lambda, row.idx = idx)
fit.bl2 <- biglasso(X       = as.big.matrix(X[idx,]),
                    y       = y[idx,],
                    family  = "cox",
                    offset  = offset[idx],
                    penalty = penalty,
                    alpha   = alpha,
                    lambda  = lambda)

as.matrix(coef(fit.gn)) - as.matrix(coef(fit.bl1))
as.matrix(coef(fit.gn)) - as.matrix(coef(fit.bl2))
as.matrix(coef(fit.bl1)) - as.matrix(coef(fit.bl2))




plot((1-fit.gn$dev.ratio)*fit.gn$nulldev ~ lambda, log = 'x', type = 'l', lwd = 2)
lines(fit.bl$loss ~ lambda, lwd = 2, col = 'red')

err <- (1-fit.gn$dev.ratio)*fit.gn$nulldev - fit.bl$loss
plot(err ~ log(lambda), lwd = 2, type = 'l')

beta <- coef(fit.gn)

idx <- 1:(n/5)
Xbigidx <- as.big.matrix(X[idx,])
dev.gn <- coxnet.deviance(x = X[idx,], y = y[idx,], offset = offset[idx], beta = beta)
dev.bl1 <- cox.deviance(X = Xbig, y = y, offset = offset, beta = beta, row.idx = idx)$dev 
dev.bl2 <- cox.deviance(X = Xbigidx, y = y[idx,], offset = offset[idx], beta = beta)$dev 

dev.bl2 - dev.bl1
dev.bl1 - dev.gn
plot(dev.gn - dev.bl ~ log(lambda), type = 'l')























