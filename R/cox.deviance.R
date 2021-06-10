#' deviance and log-likelihood calculations for fitted Cox models
#' 
#' Calculate the deviances (model, null) and log-likelihoods (saturated, model, null) 
#' associated with a set of observations and fitted coefficients.
#' 
#' @param X The design matrix, without an intercept. It must be a
#' double type \code{\link[bigmemory]{big.matrix}} object. 
#' @param y A response vector. 
#' @param beta Fitted coefficients obtained from a \code{\link[biglasso]{biglasso}} object.
#' @param row.idx Integer vector of row/observation indices to be used in computing model likelihoods.
#' @return A list with deviance and log-likelihood values calculated over the fitted coefficients.
#' \item{nulldev}{The null deviance of the model defined to be -2 * (loglik_sat - loglik_null), where the null
#' model in the context of Cox regression refers to the 0 model.} \item{dev}{Fitted deviance associated with
#' the observations and fitted coefficients. Defined to be -2 * (loglik_sat - loglik), where loglik
#' is calculated using the Breslow approximation of the Cox model likelihood.} \item{dev.ratio}{The fraction of
#' (null) deviance explained, dev.ratio = 1 - dev/nulldev. Calculated primarily as a way to quickly compare
#' outputs with those of \code{"glmnet"}.}
#' @author David Fleischer
#' 
#' Maintainer: Yaohui Zeng <yaohui.zeng@@gmail.com>
#'
#' @export cox.deviance 
cox.deviance <- function(X, y, beta, row.idx) {
  print("Hello, World!")
  y <- y[row.idx,]

  # many of the following objects are computed in the biglasso::biglasso function
  # in principle these could be passed as arguments, though I'm not sure how much of
  # a benefit that would provide (time-wise) at the cost of readability & parsimony
  tOrder      <- order(y[,1])
  d           <- as.numeric(table(y[y[,2]==1,1]))
  dtime       <- sort(unique(y[y[,2]==1,1]))
  row.idx.cox <- which(y[tOrder,1] >= min(dtime))
  d_idx       <- integer(length(row.idx.cox))
  for(i in 1:length(row.idx.cox)) d_idx[i] <- max(which(dtime <= y[tOrder[row.idx.cox[i]],1]))

  tOrig <- order(tOrder) # map from the sorted times back to the original time ordering
  # unbelievably hacky solution to omitting any censored times occurring before the first uncensored one:
  d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx)


  # # saturated likelihood (copied from glmnet)
  # ### NOTE: This may just be identical to doing -sum(d * log(d)) and is just a waste of time
  # ### to do again
  # weights <- rep(1, nrow(y)) # currently unimplemented here, but glmnet permits observation-wise weights
  # wd <- weights[y[,2] == 1]
  # tyd <- y[y[,2] == 1,1]
  # if (any(duplicated(tyd))) {
  #   wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
  # }
  # wd <- wd[wd > 0]
  # ll.sat <- -sum(wd * log(wd))
  ll.sat <- -sum(d * log(d))
  
  
  beta.T <- as(beta, "dgTMatrix")

  # Note that we are NOT calculating the R set cumulatively (hence dR for the difference/differential of
  # the R set). This allows us to save time when computing likelihoods via the fact that the cumulative
  # sets Rj are strictly decreasing as j increases. For this reason we compute the likelihoods
  # from j = length(d) ... 1 and sum the terms containing dRj as we iterate in j. See the C++ code
  # for details
  D_dR_sets <- lapply(1:length(d), function(j) {
    Dj <- d_idx2 == j
    Dj <- Dj[tOrig] & y[,2]
    dRj <- (y[,1] == y[Dj,1][1])
    list("Dj" = as.integer(which(Dj) - 1),
         "dRj" = as.integer(which(dRj) - 1))
  })

  ll.out <- loglik_cox(X@address, as.integer(row.idx - 1), beta, beta.T@i, beta.T@j, D_dR_sets, as.integer(d))
  ll <- as.numeric(ll.out)
  D <- -2 * (ll - ll.sat)

  ll.null <- 0
  # this section could be coded more elegantly, but the time it takes is trivial in
  # comparison to the model likelihood
  for (j in 1:nrow(y)) {
    if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2]
      Rj <- (y[,1] >= y[j,1])
      ll.null <- ll.null - log(sum(Rj))
    }
  }
  D0 <- -2 * (ll.null - ll.sat)

  return (list("nulldev" = D0, "dev" = D, "dev.ratio" = 1 - D/D0,
               "loglik_null" = ll.null, "loglik" = ll, "loglik_sat" = ll.sat))
  
  # ######
  # ###### IS THE ABOVE WORK EVEN NECESSARY? Can we just implement something
  # ###### simpler like below (seems to be identical to glmnet despite 
  # ###### glmnet claiming to use Breslow in the presence of ties)
  # ######
  # ### fitted likelihood
  # ll <- 0
  # for (i in 1:nrow(y)) {
  #   if (y[i,2] == 1) {
  #     j <- y[,1] >= y[i,1]
  #     
  #     ll <- ll + X[i,] %*% beta - log(colSums(exp(X[j,] %*% beta)))
  #   }
  # }
  # 
  # ### null likelihood
  # ll.null <- 0
  # for (i in 1:nrow(y)) {
  #   if (y[i,2] == 1) {
  #     j <- y[,1] >= y[i,1]
  #     
  #     # clearly this can be written more simply, but for exposition it's trivially obvious what the null model is doing
  #     #ll.null <- ll.null + X[i,] %*% rep(0, ncol(X)) - log(sum(exp(X[j,,drop=F] %*% rep(0, ncol(X))))) 
  #     # equivalent to
  #     ll.null <- ll.null - log(sum(j))
  #   }
  # }
  # 
  # ### saturated likelihood
  # weights <- rep(1, nrow(y)) # currently unimplemented here, but glmnet permits observation-wise weights
  # wd <- weights[tevent == 1]
  # tyd <- ty[tevent == 1]
  # if (any(duplicated(tyd))) {
  #   wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
  # }
  # wd <- wd[wd > 0]
  # ll.sat <- -sum(wd * log(wd))
  # 
  # D  <- -2 * (ll - ll.sat)
  # D0 <- -2 * (ll.null - ll.sat)
  # 
  # return (list("nulldev" = D0, "dev" = D, "dev.ratio" = 1 - D/D0,
  #              "loglik_null" = ll.null, "loglik" = ll, "loglik_sat" = ll.sat))
}



