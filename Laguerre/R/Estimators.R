#' Perform Laguerre Censored Quantile Regression
#'
#' @description
#' A function that estimates the model with heteroskedasticity.
#'
#' @details
#' The parameters m and m_tilde represent the order of the Laguerre series used for the positive and
#' negative parts of the error density, as described in the paper by Alexander Kreiss and Ingrid van Keilegom.
#' The parameter H represents the order of the series expansion for sigma, the heteroskedasticity function of the
#' regression model. Current options for this orthogonal expansion are Hermite, Laguerre and Legendre polynomials.
#' @param m the order of the series expansion on the positive real axis.
#' @param m_tilde the order of the series expansion on the negative real axis.
#' @param H the order of the series expansion for sigma (modulo link function).
#' @param X the design matrix. The first column must be a vector of 1's (for the intercept).
#' @param X_s the matrix of variables that sigma depends upon. Do not include an intercept.
#' @param type A vector of character values indicating the type of polynomials to be used in the expansion of sigma. Current options are "Hermite", "Laguerre", and "Legendre". If multiple variables appear in X_s, then 'type' must have as many elements as variables are in X_s. The tensor product basis will be used in this multivariate case.
#' @param Y The vector of observed values (Y = min(T,C)).
#' @param Delta The vector of censorship status (1 for observed, 0 for censored).
#' @param tau The quantile level to estimate.
#' @param starting_beta The initial values of the coefficients. Defaults to FALSE.
#' @param trials The number of initial values to try in the optimization of the log-likelihood function. Defaults to 32.
#' @param verbose Controls for output during optimization. If 0, no output is printed. If 1, a progress bar is shown indicating progress. If 2, a message is shown if an error occured during a trial. If 3, both progress bar and error alert are shown.
#' @param link The link function in sigma = link(orthogonal expansion). Options are "exp" for an exponential link, "quad" for a quadratic link, or a custom link function. A custom argument must be a list of two functions. The first entry is the link function, and the second entry is its derivative.
#' @export
laguerre_estimator_het <- function(m,m_tilde,H,X,X_s=0,type,Y,Delta,tau,starting_beta=FALSE,trials=32, verbose = 1, link="exp") {
  ## Compute Initial value for beta
  deg = H
  p <- dim(X)[2]
  H = ifelse(H==0,0,dim(Her(X_s,deg,type=type))[2] -1 )
  dims <- c(p,m,m_tilde,H)
  if(isFALSE(starting_beta)==TRUE) {
    starting_beta <- rep(0,p)
  } else if (length(starting_beta)!=p) {
    stop("Starting value has wrong dimension\n")
  }

  ##If no heteroskedasticity, use laguerre_estimator
  if(deg==0){
    est = laguerre_estimator(m,m_tilde,X,Y,Delta,tau,starting_beta=starting_beta,trials=trials)
    return(list("objective"=est$objective,"beta"=est$beta,"theta"=est$theta,"theta_tilde"=est$theta_tilde, "H"=1))
  }

  ##Compute initial estimate for Beta
  beta_no_lag=laguerre_estimator(0,0,X,Y,Delta,tau,starting_beta=starting_beta,trials=trials)$beta
  if (1==1){
    grid <- matrix(runif(trials*(H)),nrow=trials)
    H_no_lag = pi*as.matrix(grid[,1:(H)])
  }else{H_no_lag=rep(0.43,H)}


  ## If m=0 and m_tilde=0 that was it already
  if(m==0 & m_tilde==0 & deg==0) {
    return(list("objective"=-out$objective,"beta"=beta_no_lag,"theta"=1,"theta_tilde"=1, "H"=1,"X"=X,"Y"=Y,
                "X_s"=X_s,"deg"=deg, "link"=link,"type"=type,"quantile"=tau))
  }

  ## Create random grid for theta and theta_tilde
  starting_values <- matrix(0,nrow=trials,ncol=p+m+m_tilde+H)
  starting_values[,1:p] <- matrix(beta_no_lag,ncol=p,nrow=trials,byrow = TRUE)
  starting_values[,(p+m+m_tilde+1):(p+m+m_tilde+H)] = matrix(H_no_lag,ncol=H,nrow=trials,byrow = TRUE)
  ## Add trials many random points for theta and theta_tilde
  if(m!=0) {
    grid <- matrix(runif(trials*(m)),nrow=trials)
    starting_values[,(p+1):(p+m)] <- pi*as.matrix(grid[,1:(m)])
  }
  if(m_tilde!=0) {
    grid <- matrix(runif(trials*(m_tilde)),nrow=trials)
    starting_values[,(p+m+1):(p+m+m_tilde)] <- pi*as.matrix(grid[,1:(m_tilde)])
  }

  ## Do the optimization for each of the points in the grid

  result <- matrix(0,nrow=trials,ncol=1+sum(dims))
  opts <- list(algorithm="NLOPT_LD_LBFGS",print_level=0,xtol_rel=0.000001,maxeval=20000)
  if(m==0 & m_tilde==0) {
    if(verbose==1 | verbose == 3){
    cat("Optimization in progress", "\n")
    pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    }
    for(i in 1:trials) {
      out <- try(nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_beta_only2,lb=c(rep(-Inf,p),rep(0,(H))),ub=c(rep(Inf,p),rep(pi,H)),opts=opts,X=X,X_s=X_s,type=type,Y=Y,Delta=Delta,tau=tau,dims=dims, link=link, deg=deg), silent=TRUE)
      if(class(out)=="try-error"){
        result[i,1] <- Inf
        if(verbose==2|verbose==3){
          cat("Message: Trial number", i, "in optimization encountered the following error:", out,"\n")
          cat("This trial will be discarded in optimization", "\n")}
      }else{
        result[i,1] <- out$objective
        result[i,2:(sum(dims)+1)] <- out$solution
      }
      if(verbose==1 | verbose == 3){
      utils::setTxtProgressBar(pb, 100*i/trials)
        }}
    if(verbose==1|verbose==3){close(pb)}
  } else if(m==0) {
    if(verbose==1 | verbose == 3){
      cat("Optimization in progress", "\n")
      pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    }
    for(i in 1:trials) {
      out <- try(nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_beta_theta_tilde2,lb=c(rep(-Inf,p),rep(0,(m_tilde+H))),ub=c(rep(Inf,p),rep(pi,m_tilde+H)),opts=opts,X=X,X_s=X_s,type=type,Y=Y,Delta=Delta,tau=tau,dims=dims, link=link, deg=deg), silent=TRUE)
      if(class(out)=="try-error"){
        result[i,1] <- Inf
        if(verbose==2|verbose==3){
          cat("Message: Trial number", i, "in optimization encountered the following error:", out, "\n")
          cat("This trial will be discarded in optimization","\n")}
      }else{
        result[i,1] <- out$objective
        result[i,2:(sum(dims)+1)] <- out$solution
      }
      if(verbose==1 | verbose == 3){
        utils::setTxtProgressBar(pb, 100*i/trials)
      }}
    if(verbose==1|verbose==3){close(pb)}
  } else if(m_tilde==0) {
    if(verbose==1 | verbose == 3){
      cat("Optimization in progress", "\n")
      pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    }
    for(i in 1:trials) {
      out <- try(nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_beta_theta2,lb=c(rep(-Inf,p),rep(0,m+H)),ub=c(rep(Inf,p),rep(pi,m+H)),opts=opts,X=X,X_s=X_s,type=type,Y=Y,Delta=Delta,tau=tau,dims=dims, link=link,deg=deg), silent=TRUE)
      if(class(out)=="try-error"){
        result[i,1] <- Inf
        if(verbose==2|verbose==3){
          cat("Message: Trial number", i, "in optimization encountered the following error:", out,"\n")
          cat("This trial will be discarded in optimization","\n")}
      }else{
        result[i,1] <- out$objective
        result[i,2:(sum(dims)+1)] <- out$solution
      }
      if(verbose==1 | verbose == 3){
        utils::setTxtProgressBar(pb, 100*i/trials)
      }}
    if(verbose==1|verbose==3){close(pb)}
  } else {
    if(verbose==1 | verbose == 3){
      cat("Optimization in progress", "\n")
      pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
    }
    for(i in 1:trials) {
      out <- try(nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_all2,lb=c(rep(-Inf,p),rep(0,m+m_tilde+H)),ub=c(rep(Inf,p),rep(pi,m+m_tilde+H)),opts=opts,X=X,X_s=X_s,type=type,Y=Y,Delta=Delta,tau=tau,dims=dims, link=link,deg=deg), silent=TRUE)
      if(class(out)=="try-error"){
        result[i,1] <- Inf
        if(verbose==2|verbose==3){
          cat("Message: Trial number", i, "in optimization encountered the following error:", out,"\n")
          cat("This trial will be discarded in optimization","\n")}
      }else{
        result[i,1] <- out$objective
        result[i,2:(sum(dims)+1)] <- out$solution
      }
      if(verbose==1 | verbose == 3){
        utils::setTxtProgressBar(pb, 100*i/trials)
      }}
    if(verbose==1|verbose==3){close(pb)}
  }

  ## Choose the best value
  opt <- min(which(result[,1]==min(result[,1])))
  beta_est <- result[opt,2:(p+1)]
  if(m!=0) {
    theta_est <- SphericalCubature::polar2rect(1,result[opt,(p+2):(p+m+1)])
  } else {
    theta_est <- 1
  }
  if(m_tilde!=0) {
    theta_tilde_est <- SphericalCubature::polar2rect(1,result[opt,(p+m+2):(p+m+m_tilde+1)])
  } else {
    theta_tilde_est <- 1
  }

  if(H!=0) {
    H_est <- SphericalCubature::polar2rect(1,result[opt,(p+m+m_tilde+2):(p+m+m_tilde+H+1)])
  } else {
    H_est <- 1
  }
  L <- -result[opt,1]
  ## Return
  return(list("objective"=L,"beta"=beta_est,"theta"=theta_est,"theta_tilde"=theta_tilde_est, "H"=H_est,"X"=X,"Y"=Y,
              "X_s"=X_s,"deg"=deg, "link"=link,"type"=type,"quantile"=tau))
}


##########################################################################################################################
##### Laguerre estimator as in the paper. The function laguerre_estimator is called within the function "laguerre_estimator_het

## Computes the Laguerre estimator as described in the paper. For a detailed
## description of how the optimisation works we refer to the documentation
## pdf.
## Input:
##  m,m_tilde     - Integers which specify the model dimension, m=length(theta)-1
##                  and m_tilde=length(theta_tilde)-1, where theta and theta_tilde
##                  are as in likelihood
##  X             - Matrix with covariates: Every row corresponds to an observa-
##                  tion, the number of columns of X must equal length(beta), one
##                  column of X should contain only ones and correspond as such to
##                  the intercept
##  Y             - Vector of responses, length(Y) must equal nrow(X)
##  Delta         - Vector of censoring indicators, Delta[i]=1 means that observa-
##                  tion i is not censored, length(Delta)must equal length(Y)
##  tau           - Quantile of interest, element of (0,1)
##  starting_beta - If provided this is the initial value for beta used for the
##                  optimization routine, if FALSE a starting value is computed
##                  by letting m=m_tilde=0, default is ZERO
##  trials        - Default value is 32, it gives the number of random starting
##                  points for theta and theta_tilde in the optimization
## Output: List of the following four objects
##  objective   - Value of the log-likelihood at the estimated maximum
##  beta        - Estimated beta
##  theta       - Estimated theta (in Cartesian coordinates)
##  theta_tilde - Estimated theta_tilde (in Cartesian coordinates)
laguerre_estimator <- function(m,m_tilde,X,Y,Delta,tau,starting_beta=FALSE,trials=32) {
  ## Compute Initial value for beta
  p <- dim(X)[2]
  if(isFALSE(starting_beta)==TRUE) {
    starting_beta <- rep(0,p)
  } else if (length(starting_beta)!=p) {
    stop("Starting value has wrong dimension\n")
  }
  opts <- list(algorithm="NLOPT_LD_LBFGS",print_level=0,xtol_rel=0.000001,maxeval=20000)
  out <- nloptr::nloptr(x0=starting_beta,eval_f=likelihood_wrapper_beta_only,opts=opts,X=X,Y=Y,Delta=Delta,tau=tau)
  beta_no_lag <- out$solution

  ## If m=0 and m_tilde=0 that was it already
  if(m==0 & m_tilde==0) {
    return(list("objective"=-out$objective,"beta"=beta_no_lag,"theta"=1,"theta_tilde"=1))
  }

  ## Create random grid for theta and theta_tilde
  starting_values <- matrix(0,nrow=trials,ncol=p+m+m_tilde)

  ## Beta is the same for all starting values
  starting_values[,1:p] <- matrix(beta_no_lag,ncol=p,nrow=trials,byrow = TRUE)

  ## Add trials many random points for theta and theta_tilde
  if(m!=0) {
    grid <- matrix(runif(trials*m),nrow=trials)
    starting_values[,(p+1):(p+m)] <- pi*as.matrix(grid[,1:m])
  }
  if(m_tilde!=0) {
    grid <- matrix(runif(trials*m_tilde),nrow=trials)
    starting_values[,(p+m+1):(p+m+m_tilde)] <- pi*as.matrix(grid[,1:m_tilde])
  }

  ## Do the optimization for each of the points in the grid
  dims <- c(p,m,m_tilde)
  result <- matrix(0,nrow=trials,ncol=1+sum(dims))
  opts <- list(algorithm="NLOPT_LD_LBFGS",print_level=0,xtol_rel=0.000001,maxeval=20000)
  if(m==0) {
    for(i in 1:trials) {
      out <- nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_beta_theta_tilde,lb=c(rep(-Inf,p),rep(0,m_tilde)),ub=c(rep(Inf,p),rep(pi,m_tilde)),opts=opts,X=X,Y=Y,Delta=Delta,tau=tau,m_tilde=m_tilde)
      result[i,1] <- out$objective
      result[i,2:(sum(dims)+1)] <- out$solution
    }
  } else if(m_tilde==0) {
    for(i in 1:trials) {
      out <- nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_beta_theta,lb=c(rep(-Inf,p),rep(0,m)),ub=c(rep(Inf,p),rep(pi,m)),opts=opts,X=X,Y=Y,Delta=Delta,tau=tau,m=m)
      result[i,1] <- out$objective
      result[i,2:(sum(dims)+1)] <- out$solution
    }
  } else {
    for(i in 1:trials) {
      out <- nloptr::nloptr(x0=starting_values[i,],eval_f=likelihood_wrapper_all,lb=c(rep(-Inf,p),rep(0,m+m_tilde)),ub=c(rep(Inf,p),rep(pi,m+m_tilde)),opts=opts,X=X,Y=Y,Delta=Delta,tau=tau,dims=dims)
      result[i,1] <- out$objective
      result[i,2:(sum(dims)+1)] <- out$solution
    }
  }

  ## Choose the best value
  opt <- min(which(result[,1]==min(result[,1])))
  beta_est <- result[opt,2:(p+1)]
  if(m!=0) {
    theta_est <- SphericalCubature::polar2rect(1,result[opt,(p+2):(p+m+1)])
  } else {
    theta_est <- 1
  }
  if(m_tilde!=0) {
    theta_tilde_est <- SphericalCubature::polar2rect(1,result[opt,(p+m+2):(p+m+m_tilde+1)])
  } else {
    theta_tilde_est <- 1
  }
  L <- -result[opt,1]

  ## Return
  return(list("objective"=L,"beta"=beta_est,"theta"=theta_est,"theta_tilde"=theta_tilde_est))
}
