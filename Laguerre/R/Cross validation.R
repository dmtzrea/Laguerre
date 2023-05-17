#' Computes the cross-validation criterion for given model dimensions m,
#' tilde_m, and H. Note that the splitting of the data is random and happens in this
#' function. This consecutive calls of this function will yield different
#' results.
#' Input: The input variables are all identical to laguerre_cross_validation,
#'        see there for details. For this function only the model dimensions
#'        need to specified in addition.
#' Output: The value of the cross-validation criterion.
#' @export
CV_criterion <- function(m,m_tilde,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose=FALSE,link="exp") {
  ## Find random folds an save their indices
  if(nfolds==1) {
    stop("Need at least two folds, for regular estimation use laguerre_estimator")
  }
  n <- length(Y)
  K <- ceiling(n/nfolds)
  fold_indices <- matrix(0,ncol=K,nrow=nfolds)
  ind <- sample(1:n,n)
  for(i in 1:(nfolds-1)) {
    fold_indices[i,1:K] <- ind[((i-1)*K+1):(i*K)]
  }
  fold_indices[nfolds,1:(n-(nfolds-1)*K)] <- ind[((nfolds-1)*K+1):n]

  ## Compute Estimates by Leaving out the folds step by step
  beta_est <- matrix(0,nrow=nfolds,ncol=ncol(X))
  for(i in 1:nfolds) {
    ## Find indices of fold
    current_indices <- setdiff(1:n,fold_indices[i,])

    ## Compute estimate
    if(X_s !=0){
    out <- laguerre_estimator_het(m,m_tilde,H,X[current_indices,],as.matrix(X_s[current_indices]),type, Y[current_indices],Delta[current_indices],
                                  tau,starting_beta,trials,verbose=verbose, link=link)}
    else{out <- laguerre_estimator_het(m,m_tilde,H,X[current_indices,],0,type, Y[current_indices],Delta[current_indices],
                                       tau,starting_beta,trials,verbose=verbose, link=link)}
    beta_est[i,] <- out$beta
  }

  ## Compute Cross-Validation Criterion
  CVcrit <- 0
  for(i in 1:nfolds) {
    test_indices <- fold_indices[i,which(Delta[fold_indices[i,]]==1)]
    CVcrit <- CVcrit+sum(check_function(Y[test_indices]-X[test_indices,]%*%beta_est[i,],tau))/length(test_indices)
  }
  CVcrit <- CVcrit/nfolds

  return(CVcrit)
}

#' Implementation of the check function
#' @export
check_function <- function(z,tau) {
  return(z*(tau-as.numeric(z<=0)))
}


## Perform Cross-Validation as described in the paper. The function
## laguerre_estimator is repeatedly called with different values for m and
## m_tilde. The other values are always the same as specified in the call of
## laguere_cross_validation. Note that the data splitting itself happens in the
## function CV_criterion as described below. This functions is
## just an optimisation algorithm for CV_criterion.
## Input:
##  X,Y,Delta,tau - see laguerre_estimator
##  starting_beta - see laguerre_estimator
##  trials        - see laguerre_estimator
##  nfolds        - Gives the number of chunks which are build from the given
##                  data set. The dataset is thus split nfolds-1 equally sized
##                  parts and one part which is smaller than the others. These
##                  parts and then used for the cross-validation. Default value
##                  is 5.
## Output: List of three elements
##  m       - The m       yielding the best cross-validation criterion
##  m_tilde - The m_tilde yielding the best cross-validation criterion
##  est     - The output of laguerre_estimator for the best model (m,m_tilde)

#' A function that finds the best model dimensions by cross-validation. It does not perform an exhaustive search, rather, it starts
#' with the simplest model, and moves up in complexity by increasing th dimension that gives the greatest decrease in the cross-validation
#' criterion. It iterates until a local optimum is found.
#' @export
laguerre_cross_validation_het <- function(X,X_s,type,Y,Delta,tau,starting_beta=FALSE,trials=32,nfolds=5,verbose=FALSE, link="exp", h_step) {
  ## Compute initial values of CV criterion
  m <- 0
  m_tilde <- 0
  H = 0
  current_CV <- CV_criterion(m,m_tilde,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose = verbose, link)
  #left_advance_CV  <- CV_criterion(m  ,m_tilde+1,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose,link)
  #right_advance_CV <- CV_criterion(m+1,m_tilde  ,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose, link)
  #het_advance = CV_criterion(m,m_tilde, H+1, X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose, link)
  # Check if no model is already the best
  READY_FLAG <- FALSE
  #if(current_CV<=min(c(left_advance_CV,right_advance_CV,het_advance))) {
  #  READY_FLAG <- TRUE
  #}
  
  ## Increase the dimension until minimum or maximum dimension is reached
  while(!READY_FLAG) {
    ## Check in which direction to proceed in density
    left_advance_CV  <- CV_criterion(m  ,m_tilde+1,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose,link)
    right_advance_CV <- CV_criterion(m+1,m_tilde  ,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose, link)
    
    if(left_advance_CV<right_advance_CV & left_advance_CV<current_CV) {
      m_tilde <- m_tilde+1
      current_CV <- left_advance_CV
      advance_density = 0
    } else if(right_advance_CV<left_advance_CV & right_advance_CV<current_CV){
      m = m+1
      current_CV=right_advance_CV
      advance_density = 0
    } else {
      advance_density = 1
    }
    
    
    ## Check in which direction to proceed in sigma
    het_advance = CV_criterion(m,m_tilde, H+h_step, X,X_s,type,Y,Delta,tau,starting_beta,trials,nfolds,verbose, link)
    
    if(het_advance<current_CV) {
      H <- H+h_step
      current_CV <- het_advance
      advance_sigma = 0
    }else{
      advance_sigma = 1
    }
    
    
    
    if(advance_sigma*advance_density == 1 ) {
      READY_FLAG <- TRUE
    }
  }
  
  ## Compute estimate using all observations
  est <- laguerre_estimator_het(m,m_tilde,H,X,X_s,type,Y,Delta,tau,starting_beta,trials,verbose,link)
  
  return(list(m=m,m_tilde=m_tilde,H=H,est=est))
}