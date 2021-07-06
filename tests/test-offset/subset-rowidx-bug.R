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
p <- 50
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
lambda <- exp(seq(0, -5, length.out = 100))
penalty <- "enet"
alpha   <- 0.5
nfolds  <- 5
foldid  <- ceiling(sample(1:n)/n*nfolds)

idx <- which(foldid == 1)
Xbig <- as.big.matrix(X)
Xbigidx <- as.big.matrix(X[idx,])

#======================== FIT MODELS ========================#
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
fit.bl2 <- biglasso(X       = Xbigidx,
                    y       = y[idx,],
                    family  = "cox",
                    offset  = offset[idx],
                    penalty = penalty,
                    alpha   = alpha,
                    lambda  = lambda)

coef(fit.gn) - coef(fit.bl1)
coef(fit.gn) - coef(fit.bl2)
coef(fit.bl1) - coef(fit.bl2)


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























