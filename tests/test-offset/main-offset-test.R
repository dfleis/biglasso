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
p <- 10
p_nz  <- 1 # proportion of nonzerom coefs.
nb_nz <- ceiling(p_nz * p)

beta_ind <- sample(rep(c(0, 1), c(p - nb_nz, nb_nz)))
beta <- rnorm(p) * beta_ind

X <- matrix(rnorm(n * p), nrow = n)
eta <- X %*% beta
time <- sapply(eta, function(e) rexp(1, rate = exp(e)))
status <- rbinom(n, 1, 0.75)
y <- cbind("time" = time, "status" = status)

offset <- rep(0, n)
offset <- rnorm(n)*10


Xbig <- as.big.matrix(X)

# data <- coxed::sim.survdata(N = n, xvars = p, beta = beta)
# #table(data$data$y, data$data$failed)
# 
# X <- as.matrix(data$xdata)
# Xbig <- as.big.matrix(X)
# y <- cbind("time" = data$data$y, "status" = data$data$failed)
# 


#======================== FIT MODELS ========================#
### model params
penalty <- "enet"
alpha   <- 0.5
lambda  <- exp(seq(1, -7, length.out = 100))

### run models
fit.gn <- glmnet(x      = X, 
                 y      = y,
                 family = "cox", 
                 alpha  = alpha,
                 lambda = lambda)
fit.gn.off <- glmnet(x       = X, 
                     y       = y,
                     offset  = offset,
                     family  = "cox", 
                     alpha   = alpha,
                     lambda  = lambda)
fit.bl <- biglasso(X       = Xbig,
                   y       = y,
                   family  = "cox",
                   penalty = penalty,
                   alpha   = alpha,
                   lambda  = lambda)
fit.bl.off <- biglasso(X       = Xbig,
                       y       = y,
                       family  = "cox",
                       offset  = offset,
                       penalty = penalty,
                       alpha   = alpha,
                       lambda  = lambda)

beta.gn <- coef(fit.gn)
beta.gn.off <- coef(fit.gn.off)
beta.bl <- coef(fit.bl)
beta.bl.off <- coef(fit.bl.off)

sum((beta.gn - beta.bl)^2)
sum((beta.gn.off - beta.bl.off)^2)

### plots
mycol <- c(rgb(0, 0, 1, 0.75), rgb(1, 0, 0, 0.75))
mylty <- c("dashed", "dotted")

plot(NA, xlim = range(lambda), main = "Without Offset",
     ylim = range(range(beta.gn), range(beta.bl)), log = 'x', xlab = expression(log(lambda)))
grid(); abline(h = 0, lwd = 2, col = "gray50")
for (i in 1:nrow(beta.bl)) {
  lines(beta.gn[i,] ~ lambda, col = mycol[1], lty = mylty[1], lwd = 2)
  lines(beta.bl[i,] ~ lambda, col = mycol[2], lty = mylty[2], lwd = 2)
}
legend("topright", legend = c("glmnet", "biglasso"), col = mycol, lty = mylty, lwd = 2, seg.len = 2)

plot(NA, xlim = range(lambda), main = "With Offset",
     ylim = range(range(beta.gn.off), range(beta.bl.off)), log = 'x', xlab = expression(log(lambda)))
grid(); abline(h = 0, lwd = 2, col = "gray50")
for (i in 1:nrow(beta.bl.off)) {
  lines(beta.gn.off[i,] ~ lambda, col = mycol[1], lty = mylty[1], lwd = 2)
  lines(beta.bl.off[i,] ~ lambda, col = mycol[2], lty = mylty[2], lwd = 2)
}
legend("topright", legend = c("glmnet", "biglasso"), col = mycol, lty = mylty, lwd = 2, seg.len = 2)


plot(NA, xlim = range(lambda), 
     ylim = range(range(beta.gn), range(beta.bl), range(beta.gn.off), range(beta.bl.off)), 
     log = 'x', xlab = expression(log(lambda)))
grid(); abline(h = 0, lwd = 2, col = "gray50")
for (i in 1:nrow(beta.bl)) {
  lines(beta.gn[i,] ~ lambda, col = "cyan", lwd = 2, lty = 'dashed')
  lines(beta.bl[i,] ~ lambda, col = "magenta", lwd = 2, lty = 'dotted')
  lines(beta.gn.off[i,] ~ lambda, col = "red2", lwd = 2)
  lines(beta.bl.off[i,] ~ lambda, col = "blue2", lwd = 2)
}
legend("topright", legend = c("glmnet", "biglasso", "glmnet off", "biglasso off"), 
       col = c("cyan", "magenta", "red2", "blue2"), lwd = 2, 
       lty = c("dashed", "dotted", "solid", "solid"), seg.len = 2)

# err.beta <- as.matrix(beta.bl - beta.gn)
# err.beta.off <- as.matrix(beta.bl.off - beta.gn.off)
# 
# plot(NA, xlim = range(lambda), ylim = range(err.beta), log = 'x', xlab = expression(log(lambda)))
# grid(); abline(h = 0, lwd = 2, col = "gray50")
# for (i in 1:nrow(err.beta)) {
#   lines(err.beta[i,] ~ lambda, lwd = 2, col = rgb(0, 0, 0, 0.25))
# }
# plot(NA, xlim = range(lambda), ylim = range(err.beta.off), log = 'x', xlab = expression(log(lambda)))
# grid(); abline(h = 0, lwd = 2, col = "gray50")
# for (i in 1:nrow(err.beta.off)) {
#   lines(err.beta.off[i,] ~ lambda, lwd = 2, col = rgb(0, 0, 0, 0.25))
# }
# # 

