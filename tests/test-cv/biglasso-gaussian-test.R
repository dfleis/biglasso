library(biglasso)
library(glmnet)

########### genmerate data
set.seed(124)
n <- 5000
p <- 100
psparse <- 0.5
sigma <- 1

beta <- rnorm(p) * rbinom(p, 1, psparse)
X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
Xbig <- as.big.matrix(X[,-1]) # biglasso assumes the design matrix does not have the intercept column

eps <- rnorm(n, 0, sigma)
y <- X %*% beta + eps

########### model parameters
alpha  <- 0.5
lambda <- exp(seq(-2, -6, length.out = 100))
nfolds <- 5

#### raw fits
pt <- proc.time()
fit.bl <- biglasso(X = Xbig,
                   y = y,
                   family  = "gaussian",
                   penalty = "enet",
                   alpha   = alpha,
                   lambda  = lambda)
proc.time() - pt

pt <- proc.time()
fit.gn <- glmnet(x = X[,-1],
                 y = y,
                 family  = "gaussian",
                 alpha   = alpha,
                 lambda  = lambda)
proc.time() - pt

beta.bl <- as.matrix(fit.bl$beta)
beta.gn <- rbind(fit.gn$a0, as.matrix(fit.gn$beta))



#### CV
pt <- proc.time()
cvout.bl <- cv.biglasso(X = Xbig, 
                        y = y,
                        family  = "gaussian",
                        penalty = "enet",
                        alpha   = alpha,
                        lambda  = lambda,
                        nfolds  = nfolds,
                        ncores  = 1)
proc.time() - pt


pt <- proc.time()
cvout.gn <- cv.glmnet(x = X[,-1],
                      y = y, 
                      family = "gaussian",
                      alpha  = alpha,
                      lambda = lambda,
                      nfolds = nfolds,
                      foldid = cvout.bl$cv.ind,
                      parallel = F,
                      trace.it = 1)
proc.time() - pt

plot(cvout.bl)
plot(cvout.gn)

which(lambda == cvout.bl$lambda.min)
which(lambda == cvout.gn$lambda.min)







beta.bl <- as.matrix(cvout.bl$fit$beta)
beta.gn <- rbind(cvout.gn$glmnet.fit$a0, as.matrix(cvout.gn$glmnet.fit$beta))

fields::image.plot(1:p, -log(lambda), abs(beta.bl - beta.gn), col = pals::parula(128))



