#' Cross-validation for biglasso
#' 
#' Perform k-fold cross validation for penalized regression models over a grid
#' of values for the regularization parameter lambda.
#' 
#' The function calls \code{biglasso} \code{nfolds} times, each time leaving
#' out 1/\code{nfolds} of the data.  The cross-validation error is based on the
#' residual sum of squares when \code{family="gaussian"} and the binomial
#' deviance when \code{family="binomial"}.\cr \cr The S3 class object
#' \code{cv.biglasso} inherits class \code{\link[ncvreg]{cv.ncvreg}}.  So S3
#' functions such as \code{"summary", "plot"} can be directly applied to the
#' \code{cv.biglasso} object.
#' 
#' @param X The design matrix, without an intercept, as in
#' \code{\link{biglasso}}.
#' @param y The response vector, as in \code{biglasso}.
#' @param offset Optional offset vector. Only implemented for \code{family="cox"} models.
#' @param row.idx The integer vector of row indices of \code{X} that used for
#' fitting the model. as in \code{biglasso}.
#' @param family Either \code{"gaussian"}, \code{"binomial"}, \code{"cox"} or
#' \code{"mgaussian"} depending on the response. \code{"cox"} and \code{"mgaussian"}
#' are not supported yet.
#' @param eval.metric The evaluation metric for the cross-validated error and
#' for choosing optimal \code{lambda}. "default" for linear regression is MSE
#' (mean squared error), for logistic regression is misclassification error.
#' "MAPE", for linear regression only, is the Mean Absolute Percentage Error.
#' @param grouped Currently is only applicable to Cox models, otherwise ignored.
#' Mirrors the \code{grouped} parameter used by \code{glmnet}. If \code{grouped}=TRUE 
#' then the partial likelihood for the K-th fold is calculated by subtracting
#' the (log) partial likelihood evaluated on the full data set from the 
#' partial likelihood evaluated on (K-1)/K data set. If \code{grouped}=FALSE
#' then the partial likelihood is computed only on the K-th fold.
#' @param ncores The number of cores to use for parallel execution across a
#' cluster created by the \code{parallel} package. (This is different from
#' \code{ncores} in \code{\link{biglasso}}, which is the number of OpenMP
#' threads.)
#' @param ... Additional arguments to \code{biglasso}.
#' @param nfolds The number of cross-validation folds.  Default is 5.
#' @param seed The seed of the random number generator in order to obtain
#' reproducible results.
#' @param cv.ind Which fold each observation belongs to.  By default the
#' observations are randomly assigned by \code{cv.biglasso}.
#' @param trace If set to TRUE, cv.biglasso will inform the user of its
#' progress by announcing the beginning of each CV fold.  Default is FALSE.
#' @return An object with S3 class \code{"cv.biglasso"} which inherits from
#' class \code{"cv.ncvreg"}.  The following variables are contained in the
#' class (adopted from \code{\link[ncvreg]{cv.ncvreg}}).  \item{cve}{The error
#' for each value of \code{lambda}, averaged across the cross-validation
#' folds.} \item{cvse}{The estimated standard error associated with each value
#' of for \code{cve}.} \item{lambda}{The sequence of regularization parameter
#' values along which the cross-validation error was calculated.}
#' \item{fit}{The fitted \code{biglasso} object for the whole data.}
#' \item{min}{The index of \code{lambda} corresponding to \code{lambda.min}.}
#' \item{lambda.min}{The value of \code{lambda} with the minimum
#' cross-validation error.} \item{null.dev}{The deviance for the intercept-only
#' model.} \item{pe}{If \code{family="binomial"}, the cross-validation
#' prediction error for each value of \code{lambda}.} \item{cv.ind}{Same as
#' above.}
#' @author Yaohui Zeng and Patrick Breheny
#' 
#' Maintainer: Yaohui Zeng <yaohui.zeng@@gmail.com>
#' @seealso \code{\link{biglasso}}, \code{\link{plot.cv.biglasso}},
#' \code{\link{summary.cv.biglasso}}, \code{\link{setupX}}
#' @examples
#' \dontrun{
#' ## cv.biglasso
#' data(colon)
#' X <- colon$X
#' y <- colon$y
#' X.bm <- as.big.matrix(X)
#' 
#' ## logistic regression
#' cvfit <- cv.biglasso(X.bm, y, family = 'binomial', seed = 1234, ncores = 2)
#' par(mfrow = c(2, 2))
#' plot(cvfit, type = 'all')
#' summary(cvfit)
#' }
#' 
#' @export cv.biglasso
#' 
cv.biglasso <- function(X, y, offset = NULL, row.idx = 1:nrow(X), 
                        family = c("gaussian", "binomial", "cox", "mgaussian"),
                        eval.metric = c("default", "MAPE"),
                        grouped = TRUE,
                        ncores = parallel::detectCores(), ...,
                        nfolds = 5, seed, cv.ind, trace = FALSE) {
  
  family <- match.arg(family)
  if(!family %in% c("gaussian", "binomial", "cox")) stop("CV method for this family not supported yet.")
  
  #TODO: 
  #   system-specific parallel: Windows parLapply; others: mclapply
  eval.metric <- match.arg(eval.metric)

  max.cores <- parallel::detectCores()
  if (ncores > max.cores) {
    cat("The number of cores specified (", ncores, ") is larger than the number of avaiable cores (", max.cores, "). Use ", max.cores, " cores instead! \n", sep = "")
    ncores = max.cores
  }
  
  fit <- biglasso(X = X, y = y, offset = offset, row.idx = row.idx, family = family, ncores = ncores, ...)
  n <- fit$n
  if (fit$family == "cox") {
    E <- matrix(NA, nrow=nfolds, ncol=length(fit$lambda))
    w <- rep(NA, nfolds) # CV weights (w[i] = # of uncensored obs. in the i-th test set, copied from glmnet)
  } else {
    E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))
    # y <- fit$y # this would cause error if eval.metric == "MAPE"
  }
  
  if (fit$family == 'binomial') {
    PE <- E
  }

  if (!missing(seed)) set.seed(seed)
  if (missing(cv.ind)) {
    if (fit$family=="binomial" & (min(table(y)) > nfolds)) {
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
      cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
      cv.ind <- numeric(n)
      cv.ind[y==1] <- cv.ind1
      cv.ind[y==0] <- cv.ind0
    } else {
      cv.ind <- ceiling(sample(1:n)/n*nfolds)
    }
  }
  
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$family <- fit$family # I think omitting this may have been a bug?
  
  parallel <- FALSE
  if (ncores > 1) {
    cluster <- parallel::makeCluster(ncores)
    if (!("cluster" %in% class(cluster))) stop("cluster is not of class 'cluster'; see ?makeCluster")
    parallel <- TRUE
    ## pass the descriptor info to each cluster ##
    xdesc <- bigmemory::describe(X)
    parallel::clusterExport(cluster, c("cv.ind", "xdesc", "y", "cv.args", 
                                       "parallel", "eval.metric"), 
                            envir=environment())
    parallel::clusterCall(cluster, function() {
     
      require(biglasso)
      # require(bigmemory)
      # require(Matrix)
      # dyn.load("~/GitHub/biglasso.Rcheck/biglasso/libs/biglasso.so")
      # source("~/GitHub/biglasso/R/biglasso.R")
      # source("~/GitHub/biglasso/R/predict.R")
      # source("~/GitHub/biglasso/R/loss.R")
    })
    fold.results <- parallel::parLapply(cl = cluster, X = 1:nfolds, fun = cvf, XX = xdesc, 
                                        y = y, offset = offset, eval.metric = eval.metric, 
                                        cv.ind = cv.ind, cv.args = cv.args, 
                                        grouped = grouped,
                                        parallel = parallel)
    parallel::stopCluster(cluster)
  }

  for (i in 1:nfolds) {
    if (parallel) {
      res <- fold.results[[i]]
    } else {
      if (trace) cat("Starting CV fold #", i, sep="", "\n")
      res <- cvf(i, X, y, offset, eval.metric, cv.ind, cv.args, grouped)
    }
    
    if (fit$family=="cox") {
      E[i, 1:res$nl] <- res$loss # 1:res$nl to account for saturated values that had been removed
      w[i]           <- res$weight
    } else {
      E[cv.ind == i, 1:res$nl] <- res$loss
      if (fit$family == "binomial") PE[cv.ind == i, 1:res$nl] <- res$pe
      Y[cv.ind == i, 1:res$nl] <- res$yhat
    }
  }

  ## Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  if (fit$family != "cox") Y <- Y[,ind]
  lambda <- fit$lambda[ind]
  
  ## Return
  if (fit$family == "cox") {
    cve <- apply(E, 2, weighted.mean, w = w)
    cvse <- sqrt(apply(scale(E, cve, FALSE)^2, 2, weighted.mean, w = w)/(nfolds - 1))
  } else {
    cve <- apply(E, 2, mean)
    cvse <- apply(E, 2, sd) / sqrt(n)
  }
  min <- which.min(cve)

  ## include values for the 1-se rule
  cve.1se   <- (cve + cvse)[min]
  lambda.1se <- max(lambda[cve <= cve.1se], na.rm = T)
  
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, min=min, lambda.min=lambda[min],
              cve.1se = cve.1se, lambda.1se = lambda.1se,
              cv.ind = cv.ind,
              eval.metric = eval.metric, family=fit$family)
  if (fit$family=="cox") {
    ## NOTE: Do I need to calculate the entire deviance again? Can I just pass some of the 
    # results through from the cvf() call or calculate from the first call of biglasso?
    val$null.dev <- cox.deviance(X = X, y = y, offset = offset, beta = fit$beta, row.idx = 1:nrow(y))$nulldev
  } else {
    val$null.dev <- mean(loss.biglasso(y, rep(mean(y), n), 
                                       fit$family, eval.metric = eval.metric))
  }
  
  if (fit$family=="binomial") {
    pe <- apply(PE, 2, mean)
    val$pe <- pe[is.finite(pe)]
  }
  # if (returnY) val$Y <- Y
  structure(val, class=c("cv.biglasso", "cv.ncvreg"))
}

cvf <- function(i, XX, y, offset, eval.metric, cv.ind, cv.args, grouped, parallel=FALSE) {
  # reference to the big.matrix by descriptor info
  if (parallel) {
    XX <- attach.big.matrix(XX)
  }
  
  cv.args$X <- XX
  cv.args$y <- y
  cv.args$offset <- offset
  cv.args$row.idx <- which(cv.ind != i)
  cv.args$warn <- FALSE
  cv.args$ncores <- 1

  idx.test <- which(cv.ind == i)
  fit.i <- do.call("biglasso", cv.args)

  if (fit.i$family=="cox") {
    ##
    ## TODO: check if size of fold makes sense (copy glmnet's behaviour)
    ## Warning message:
    ## Option grouped=TRUE enforced for cv.coxnet, since < 10 observations per fold 
    ##
    ## It's probably also worth creating an error/warning when too many censored obs.
    ## are in a single fold (I believe glmnet does this as well)
    ##
    wt <- sum(y[idx.test, 2]) # NOTE: If all obs. are censored in a test set then this will lead to a div by zero
    if (grouped) { # "V&VH cross-validation error" (default setting in glmnet)
      plfull   <- cox.deviance(X = XX, y = y, offset = offset, beta = fit.i$beta, row.idx = 1:nrow(y))
      plminusk <- cox.deviance(X = XX, y = y, offset = offset, beta = fit.i$beta, row.idx = which(cv.ind != i))
      loss <- (plfull$dev - plminusk$dev)/wt
    } else { # "basic cross-validation error" (alternative setting in glmnet)
      plk <- cox.deviance(X = XX, y = y, offset = offset, beta = fit.i$beta, row.idx = idx.test)
      loss <- plk$dev/wt
    }
    
  } else {
    y2 <- y[cv.ind==i]
    yhat <- matrix(predict(fit.i, XX, row.idx = idx.test, type="response"), length(y2))
    loss <- loss.biglasso(y2, yhat, fit.i$family, eval.metric = eval.metric)
  }
  
  pe <- if (fit.i$family=="binomial") {(yhat < 0.5) == y2} else NULL
  out <- list(loss=loss, pe=pe, nl=length(fit.i$lambda))
  if (fit.i$family == "cox") {
    out$weight <- wt
  } else {
    out$yhat <- yhat
  }
  return (out)
}

