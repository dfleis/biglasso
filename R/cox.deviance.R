#' deviance and log-likelihood calculations for fitted Cox models
#' 
#' Calculate the deviances (model, null) and log-likelihoods (saturated, model, null) 
#' associated with a set of observations and fitted coefficients.
#' 
#' @param X The design matrix, without an intercept. It must be a
#' double type \code{\link[bigmemory]{big.matrix}} object. 
#' @param y A response vector. 
#' @param offset Optional offset vector.
#' @param beta Fitted coefficients obtained from a \code{\link[biglasso]{biglasso}} object. Expected to be 
#' a matrix (or something coercible to a \code{"dgTMatrix"} object) with ncols corresponding to the number
#' of lambdas used in the penalization path and nrows corresponding to the number of model covariates.
#' @param row.idx Integer vector of row/observation indices to be used in computing model likelihoods.
#' Default set to \code{1:nrow(X)}.
#' @return A list with deviance and log-likelihood values calculated over the fitted coefficients.
#' \item{nulldev}{The null deviance of the model defined to be -2 * (loglike_sat - loglike_null), where the null
#' model in the context of Cox regression refers to the 0 model.} \item{dev}{Fitted deviance associated with
#' the observations and fitted coefficients. Defined to be -2 * (loglike_sat - loglike), where loglike
#' is calculated using the Breslow approximation of the Cox model likelihood.} \item{dev.ratio}{The fraction of
#' (null) deviance explained, dev.ratio = 1 - dev/nulldev. Calculated primarily as a way to quickly compare
#' outputs with those of \code{"glmnet"}.}
#' @author David Fleischer
#' 
#' Maintainer: Yaohui Zeng <yaohui.zeng@@gmail.com>
#'
#' @export cox.deviance 
cox.deviance <- function(X, y, offset = NULL, beta, row.idx = 1:nrow(X)) {
  y <- y[row.idx,]
  
  if (is.null(offset)) offset <- as.double(rep(0, nrow(X)))
  
  d     <- as.numeric(table(y[y[,2]==1,1]))
  dtime <- sort(unique(y[y[,2]==1,1]))
  
  ### saturated likelihood
  # NOTE: If we wish to include weighted observations, check the saturated likelihood calculations
  # from glmnet::coxnet.deviance for an elegant solution
  # https://github.com/cran/glmnet/blob/f4fc95ab49efaad9b6e1728a7c840bc6159501dc/R/coxnet.deviance.R
  ll.sat <- -sum(d * log(d))
  
  ### null likelihood (from the 0 model)
  # this could be much more efficient/elegant, but the runtime is trivial comparison to the model likelihood
  # to be honest, we probably could just call the C++ function loglik_cox with beta set to 0, but I haven't
  # yet verified it works
  ll.null <- 0
  for (j in 1:length(d)) {
    Dset <- (y[,1] == dtime[j]) & (y[,2] == 1)
    term1 <- sum(offset[Dset])
    
    Rset <- (y[,1] >= dtime[j]) 
    term2 <- d[j] * log(sum(exp(offset[Rset])))
    
    ll.null <- ll.null + term1 - term2
  }
  D0 <- -2 * (ll.null - ll.sat) 

  ### fitted likelihoods
  tOrder      <- order(y[,1])
  row.idx.cox <- which(y[tOrder,1] >= min(dtime))
  d_idx       <- integer(length(row.idx.cox))
  for(i in 1:length(row.idx.cox)) d_idx[i] <- max(which(dtime <= y[tOrder[row.idx.cox[i]],1]))

  tOrig <- order(tOrder) # map from the sorted times back to the original time ordering
  # unbelievably hacky solution to omitting any censored times occurring before the first uncensored one:
  d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx)

  beta.T <- as(beta, "dgTMatrix") # will break if beta is a vector and not a matrix
  
  # Note that we are NOT calculating the R set cumulatively (hence dR for the difference/differential of
  # the R set). This allows us to save time when computing likelihoods via the fact that the cumulative
  # sets Rj are strictly decreasing as j increases. For this reason we compute the likelihoods
  # backwards, from j = length(d) ... 1, and sum the terms containing dRj as we iterate in j.
  # See the C++ code for details
  D_dR_sets <- lapply(1:length(d), function(j) {
    dRset <- (d_idx2 == j)[tOrig]
    Dset  <- dRset & y[,2]
    list("D"  = as.integer(which(Dset) - 1), 
         "dR" = as.integer(which(dRset) - 1))
    })
  
  #return (list("nulldev" = D0, "loglike_null" = ll.null, "loglike_sat" = ll.sat))  
  ll.out <- loglik_cox(X@address, 
                       offset, 
                       as.integer(row.idx - 1),
                       beta,
                       beta.T@i, 
                       beta.T@j,
                       D_dR_sets, 
                       as.integer(d))
  ll <- as.numeric(ll.out)
  D  <- -2 * (ll - ll.sat)
  return (list("nulldev" = D0, "dev" = D, "dev.ratio" = 1 - D/D0,
               "loglike_null" = ll.null, "loglike" = ll, "loglike_sat" = ll.sat))
}

# cox.deviance <- function(X, y, offset = NULL, beta, row.idx = 1:nrow(X)) {
#   y <- y[row.idx,]
#   if (is.null(offset)) offset <- as.double(rep(0, nrow(X)))
# 
#   tOrder      <- order(y[,1])
#   d           <- as.numeric(table(y[y[,2]==1,1]))
#   dtime       <- sort(unique(y[y[,2]==1,1]))
#   row.idx.cox <- which(y[tOrder,1] >= min(dtime))
#   d_idx       <- integer(length(row.idx.cox))
#   for(i in 1:length(row.idx.cox)) d_idx[i] <- max(which(dtime <= y[tOrder[row.idx.cox[i]],1]))
# 
#   tOrig <- order(tOrder) # map from the sorted times back to the original time ordering
#   # unbelievably hacky solution to omitting any censored times occurring before the first uncensored one:
#   d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx)
#   
#   beta.T <- as(beta, "dgTMatrix") # will break if beta is a vector and not a matrix
#   
#   ### saturated likelihood 
#   ll.sat <- -sum(d * log(d))
#   
#   ### fitted likelihood
#   # Note that we are NOT calculating the R set cumulatively (hence dR for the difference/differential of
#   # the R set). This allows us to save time when computing likelihoods via the fact that the cumulative
#   # sets Rj are strictly decreasing as j increases. For this reason we compute the likelihoods
#   # backwards, from j = length(d) ... 1, and sum the terms containing dRj as we iterate in j. 
#   # See the C++ code for details
#   D_dR_sets <- lapply(1:length(d), function(j) {
#     dRj <- (d_idx2 == j)[tOrig]
#     Dj  <- dRj & y[,2]
#     list("Dj" = as.integer(which(Dj) - 1),
#          "dRj" = as.integer(which(dRj) - 1))
#   })
#   
#   ll.out <- loglik_cox(X@address, offset, as.integer(row.idx - 1), 
#                        beta, beta.T@i, beta.T@j, D_dR_sets, as.integer(d))
#   ll <- as.numeric(ll.out)
#   D <- -2 * (ll - ll.sat)
# 
#   ### null likelihood (from the 0 model)
#   ll.null <- 0
#   # this section could be coded more elegantly, but the time it takes is trivial in
#   # comparison to the model likelihood
#   for (j in 1:nrow(y)) {
#     if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2]
#       Rj <- (y[,1] >= y[j,1])
#       ll.null <- ll.null + offset[j] - log(sum(exp(offset[Rj])))
#     }
#   }
#   
#   D0 <- -2 * (ll.null - ll.sat)
# 
#   return (list("nulldev" = D0, "dev" = D, "dev.ratio" = 1 - D/D0,
#                "loglike_null" = ll.null, "loglike" = ll, "loglike_sat" = ll.sat))
# 
#   # ###### IS THE ABOVE WORK EVEN NECESSARY? 
#   # ###### Can we just implement something simpler (see like below 
#   # ###### -- seems to be identical to glmnet despite glmnet stating
#   # ###### it uses Breslow's approximation in the presence of ties)
#   # ######
#   # ###### It might be a good idea to implement it in a separate function
#   # ###### and decide later?
#   # ###### Obviously this would have to have a C++ section as well to
#   # ###### account for big.matrix objects.
#   #
#   # ### fitted likelihood
#   # ll <- 0
#   # for (i in 1:nrow(y)) {
#   #   if (y[i,2] == 1) {
#   #     j <- y[,1] >= y[i,1]
#   #     
#   #     ll <- ll + X[i,] %*% beta - log(colSums(exp(X[j,] %*% beta)))
#   #   }
#   # }
#   # 
#   # ### null likelihood
#   # ll.null <- 0
#   # for (i in 1:nrow(y)) {
#   #   if (y[i,2] == 1) {
#   #     j <- y[,1] >= y[i,1]
#   #     
#   #     # clearly this can be written more simply, but for exposition it's trivially obvious what the null model is doing
#   #     #ll.null <- ll.null + X[i,] %*% rep(0, ncol(X)) - log(sum(exp(X[j,,drop=F] %*% rep(0, ncol(X))))) 
#   #     # equivalent to
#   #     ll.null <- ll.null - log(sum(j))
#   #   }
#   # }
#   # 
#   # ### saturated likelihood
#   # # this version is from glmnet::coxnet.deviance(). I believe it's identical to -sum(d * log(d))
#   # # when weights are = 1, but clearly allows for more flexible modelling when given user-supplied weights
#   # weights <- rep(1, nrow(y)) # currently unimplemented here, but glmnet permits observation-wise weights
#   # wd <- weights[tevent == 1]
#   # tyd <- ty[tevent == 1]
#   # if (any(duplicated(tyd))) {
#   #   wd <- tapply(wd, tyd, sum) # counts of unique failure times for uncens. obs (equiv. to table(y[y[,2]==1,1]))
#   # }
#   # wd <- wd[wd > 0]
#   # ll.sat <- -sum(wd * log(wd))
#   # 
#   # D  <- -2 * (ll - ll.sat)
#   # D0 <- -2 * (ll.null - ll.sat)
#   # 
#   # return (list("nulldev" = D0, "dev" = D, "dev.ratio" = 1 - D/D0,
#   #              "loglik_null" = ll.null, "loglik" = ll, "loglik_sat" = ll.sat))
# }



