library(biglasso)
library(Rcpp)


set.seed(124)
### GENERATE DATA 
no <- 1000
nx <- 5
dat <- coxed::sim.survdata(N = no, xvars = nx)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

penalty <- "enet"
alpha   <- 0.5
lambda  <- exp(seq(-1, -6, length.out = 100)) #c(0.1, 0.01)

fit.bl <- biglasso::biglasso(X = Xbig,
                             y = y,
                             family = "cox",
                             penalty = penalty,
                             alpha   = alpha,
                             lambda  = lambda)
fit.gn <- glmnet::glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)
beta <- fit.bl$beta
lambda <- fit.bl$lambda

pt <- proc.time()
out <- biglasso::deviance.cox(Xbig, y, beta, 1:nrow(y))
proc.time() - pt
out$nulldev
fit.gn$nulldev

d1 <- out$dev
d2 <- fit.bl$loss
d3 <- (1-fit.gn$dev.ratio)*fit.gn$nulldev

plot(d1 ~ lambda, log = 'x', type = 'l', lwd = 3)
lines(d2 ~ lambda, col = 'red', lwd = 3)
lines(d3 ~ lambda, col = 'blue', lwd = 3)

plot(abs(d1 - d3) ~ lambda, log = 'x', type = 'l', lwd = 3)



m <- rbind(c(out$dev, out$nulldev), 
           c((1-fit.gn$dev.ratio)*fit.gn$nulldev, fit.gn$nulldev))
rownames(m) <- c("test", "glmnet")
colnames(m) <- c("D", "D0")
m


eta1 <- as.matrix(out)
eta2 <- X[rowidx,] %*% as.matrix(beta)
eta1-eta2

colSums(X[1:5,] %*% as.matrix(beta))
colSums(X[1:5,]) %*% as.matrix(beta)


#RcppExport SEXP get_eta(SEXP xP, SEXP row_idx_, SEXP beta, SEXP idx_p, SEXP idx_l) {
  
#cv.bl <- biglasso::cv.biglasso(X = Xbig, y = y, penalty = "enet", alpha = 0.5, nfolds = 3, trace = T, ncores = 1)


# res <- .Call("cdfit_gaussian_ada_edpp_ssr", X@address, yy, as.integer(row.idx-1),
#              lambda, as.integer(nlambda), as.integer(lambda.log.scale),
#              lambda.min, alpha,
#              as.integer(user.lambda | any(penalty.factor==0)),
#              eps, as.integer(max.iter), penalty.factor,
#              as.integer(dfmax), as.integer(ncores), update.thresh, as.integer(verbose),
#              PACKAGE = 'biglasso')