#' Internal biglasso functions
#' 
#' Internal biglasso functions
#' 
#' These are not intended for use by users.\code{loss.biglasso} calculates the 
#' value of the loss function for the given predictions (used for cross-validation).
#' 
#' @param X The design matrix, without an intercept. It must be a
#' double type \code{\link[bigmemory]{big.matrix}} object. 
#' @param y A response vector. 
#' @param beta Fitted coefficients obtained from a \code{\link[biglasso]{biglasso}} object.
#' @param row.idx Integer vector of row/observation indices that should be used to compute model likelihoods.
#' @return todo...
#' @author David Fleischer
#' 
#' Maintainer: Yaohui Zeng <yaohui.zeng@@gmail.com>
#'
#' @export deviance.cox 
deviance.cox <- function(X, y, beta, row.idx) {
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
  # unbelievably hacky solution to omitting any censored times occurring before the first uncensored one
  d_idx2 <- c(rep(-999, length(tOrder) - length(row.idx.cox)), d_idx)
  
  ll.sat <- -sum(d * log(d))
  
  beta.T <- as(beta, "dgTMatrix") 
  #eta <- get_eta(X@address, as.integer(row.idx-1), beta, beta.T@i, beta.T@j)
  
  ##### TO DO: might be more efficient to wrap this in some sort of *apply function
  # or put the entire thing into C++, but for a first pass it'll be a loop in R
  #ll <- rep(0, beta@Dim[2])
  # print("ll loop")
  # for (j in 1:length(d)) {
  #   Dj <- d_idx2 == j # indices of the death/failures corresponding to the j-th unique failure time
  #   Dj <- Dj[tOrig] & y[,2] # hacky solution to remove censored obs & place the indices in the original ordering
  # 
  #   XDbeta_raw <- get_eta(X@address, as.integer(which(Dj) - 1), beta, beta.T@i, beta.T@j)
  #   XDbeta <- as.matrix(XDbeta_raw)
  # 
  #   term1 <- colSums(XDbeta)
  # 
  #   Rj <- (y[,1] >= y[Dj,1][1]) # all times in Dj will be the same, so we pick the first
  #   ############ THIS IS THE SLOW STEP:
  #   XRbeta_raw <- get_eta(X@address, as.integer(which(Rj) - 1), beta, beta.T@i, beta.T@j)
  # 
  #   XRbeta <- as.matrix(XRbeta_raw)
  #   term2 <- d[j] * log(colSums(exp(XRbeta)))
  # 
  #   ll <- ll + term1 - term2
  # }
  # D <- -2 * (ll - ll.sat)
  D_R_sets <- lapply(1:length(d), function(j) {
    Dj <- d_idx2 == j
    Dj <- Dj[tOrig] & y[,2]
    Rj <- (y[,1] >= y[Dj,1][1])
    list("Dj" = Dj, "Rj" = Rj, "dj" = d[j])
  })
  
  ll.list <- sapply(D_R_sets, function(DRs) {
    #term1 <- colSums(X[DRs$Dj,,drop=F]) %*% beta
    #term2 <- DRs$d * log(colSums(exp(X[DRs$Rj,,drop=F] %*% beta)))
    XDbeta_raw <- get_eta(X@address, as.integer(which(DRs$Dj) - 1), beta, beta.T@i, beta.T@j)
    XDbeta <- as.matrix(XDbeta_raw)
    term1 <- colSums(XDbeta)
    
    XRbeta_raw <- get_eta(X@address, as.integer(which(DRs$Rj) - 1), beta, beta.T@i, beta.T@j)
    XRbeta <- as.matrix(XRbeta_raw)
    term2 <- DRs$dj * log(colSums(exp(XRbeta)))

    term1 - term2
  })
  D <- -2 * (rowSums(ll.list) - ll.sat)
  
  
   
  ll.null <- 0
  for (j in 1:nrow(y)) {
    if (y[j,2] == 1) { # equiv to multiplying the final difference (term1 - term2) by y[j,2]
      Rj <- (y[,1] >= y[j,1])
      ll.null <- ll.null - log(sum(Rj))
    }
  }
  D0 <- -2 * (ll.null - ll.sat)
  
  list("nulldev" = D0, "dev" = D)#, "loglik_null" = ll.null, "loglik" = ll, "loglik_sat" = ll.sat)
}



