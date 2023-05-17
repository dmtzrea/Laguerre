likelihood2 <- function(tau,beta,theta,theta_tilde,H,X,X_s,type,Y,Delta,derivative=FALSE, link="exp",deg) {

  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  he = length(H)
  ### Computing Hermite polynomials
  ### These are the terms in each variable up to the specified degree. If only one
  ### variable is considered for heteroskedasticity, then no cross-products will be present.
  Her = Her(X_s,deg, type=type)
if(length(link) == 1){
  if (link == "exp"){
    sigma = exp(Her%*%H)/exp(1)
    dsigma = sigma
  }

  if (link == "quad"){
    sigma = (Her%*%H)^2
    dsigma = 2*(Her%*%H)
  }
}else{
    sigma = as.vector(unlist(lapply(link, function(f) f(Her%*%H))[1]))
    dsigma = as.vector(unlist(lapply(link, function(f) f(Her%*%H))[2]))
  }

  z <- (Y-X%*%as.matrix(beta))/sigma
  n <- length(z)

  ## Split data in censored and uncensored part and postive and negative
  up_index <- (Delta==1) & (z>=0)
  um_index <- (Delta==1) & (z <0)
  cp_index <- (Delta==0) & (z>=0)
  cm_index <- (Delta==0) & (z <0)

  z_uncensored <- z[Delta==1]
  z_censored   <- z[Delta==0]

  zup <- z_uncensored[z_uncensored>=0]

  zum <- z_uncensored[z_uncensored <0]

  zcp <- z_censored[z_censored>=0]
  zcm <- z_censored[z_censored <0]

  s_censored = sigma[Delta==0]
  s_uncensored = sigma[Delta==1]

  scp = s_censored[z_censored>=0]
  scn = s_censored[z_censored<0]
  sup = s_uncensored[z_uncensored>=0]
  sun = s_uncensored[z_uncensored<0]

  ## Compute the matrix of binomial coefficients
  B       <- matrix(0,ncol=m+1      ,nrow=m+1)
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)

  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }

  #### Compute the density for the uncensored observations
  ## Compute the powers of z
  Zp <- matrix(1,ncol=m+1,nrow=length(zup))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zp[,i] <- Zp[,i-1]*(-tau*zup)/(i-1)
    }
  }
  Zm <- matrix(1,ncol=m_tilde+1,nrow=length(zum))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zm[,i] <- Zm[,i-1]*((1-tau)*zum)/(i-1)
    }
  }

  ## Compute Laguerre Polynomials
  Lp <- Zp%*%t(B)
  Lm <- Zm%*%t(B_tilde)

  ## Compute the density
  f <- rep(0,length(z_uncensored))
  f[z_uncensored>=0] <- (1/sup)*(1-tau)*tau*exp(-tau*zup)*(Lp%*%theta)^2
  f[z_uncensored< 0] <- (1/sun)*tau*(1-tau)*exp((1-tau)*zum)*(Lm%*%theta_tilde)^2

  #### Compute the distribution function for censored observations
  ## Compute exponential integrals
  Iupper <- matrix(0,ncol=2*m_tilde+1,nrow=length(zcm))
  Ilower <- matrix(0,ncol=2*m      +1,nrow=length(zcp))

  Iupper[,1] <-   exp((1-tau)*zcm)
  Ilower[,1] <- 1-exp(   -tau*zcp)

  if(m_tilde>=1) {
    for(k in 2:(2*m_tilde+1)) {
      Iupper[,k] <- (-(1-tau)*zcm)^(k-1)*exp((1-tau)*zcm)+(k-1)*Iupper[,k-1]
    }
  }
  if(m>=1) {
    for(k in 2:(2*m+1)) {
      Ilower[,k] <- -(tau*zcp)^(k-1)*exp(-tau*zcp)+(k-1)*Ilower[,k-1]
    }
  }

  ## Compute help vectors
  v       <- (-1)^(0:m)      /factorial(0:m)
  v_tilde <- (-1)^(0:m_tilde)/factorial(0:m_tilde)
  h       <- (t(theta)%*%B            )*v
  h_tilde <- (t(theta_tilde)%*%B_tilde)*v_tilde

  ## Compute distribution function
  F <- rep(0,length(z_censored))
  neg_index <- z_censored<0
  pos_index <- z_censored>=0
  F[pos_index] <- tau
  for(i1 in 0:m_tilde) {
    for(i2 in 0:m_tilde) {
      F[neg_index] <- F[neg_index]+tau*h_tilde[i1+1]*h_tilde[i2+1]*Iupper[,i1+i2+1]
    }
  }
  for(i1 in 0:m) {
    for(i2 in 0:m) {
      F[pos_index] <- F[pos_index]+(1-tau)*h[i1+1]*h[i2+1]*Ilower[,i1+i2+1]
    }
  }

  ## Compute log-Likelihood
  L <- sum(log(f))+sum(log(1-F))

  #### Compute the derivative or return
  if(derivative==FALSE) {
    return(L/n)
  }

  ## Derivative of log-density with respect to sigma (H)

  if(sum(up_index)>=1) {
    if(m==0) {
      logfgradp <- Her[up_index,]*(((zup*tau/sigma[up_index])-(1/sigma[up_index]))*dsigma[up_index])
    } else {
      if(m==1) {
        helper <-  t(t(Zp[,1:m      ]))
      } else if(dim(Zp)[1]==1) {
        helper <-  t(Zp[,1:m      ])
      } else {
        helper <-  Zp[,1:m      ]
      }
      logfgradp <-Her[up_index,]*((tau*(zup/sigma[up_index])*as.numeric(1+2*helper%*%(t(B      )%*%theta      )
                                                                        [2:(m      +1)]/(Lp%*%theta      )) -
                                     (1/sigma[up_index]))*dsigma[up_index])
    }
    if(sum(up_index)==1) {
      logfgradp <- matrix(logfgradp,nrow=1)
    }
  } else {
    logfgradp <- matrix(0,nrow=1,ncol=length(H))
  }

  if(sum(um_index)>=1) {
    if(m_tilde==0) {
      logfgradm <- Her[um_index,]*(((zum*tau/sigma[um_index])-(1/sigma[um_index]))*dsigma[um_index])
    } else {
      if(m_tilde==1) {
        helper <- t(t(Zm[,1:m_tilde]))
      } else if(dim(Zm)[1]==1) {
        helper <- t(Zm[,1:m_tilde])
      } else {
        helper <- Zm[,1:m_tilde]
      }
      logfgradm <- Her[um_index,]*((-(1-tau)*(zum/sigma[um_index])*as.numeric(1+2*helper%*%(t(B_tilde)%*%theta_tilde)
                                                                              [2:(m_tilde+1)]/(Lm%*%theta_tilde)) -
                                      (1/sigma[um_index]))*dsigma[um_index])
    }
    if(sum(um_index)==1) {
      logfgradm <- matrix(logfgradm,nrow=1)
    }
  } else {
    logfgradm <- matrix(0,nrow=1,ncol=length(H))
  }

  ## Compute the density for the censored observations
  ## Compute the powers of z
  Zpc <- matrix(1,ncol=m+1,nrow=length(zcp))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zpc[,i] <- Zpc[,i-1]*(-tau*zcp)/(i-1)
    }
  }
  Zmc <- matrix(1,ncol=m_tilde+1,nrow=length(zcm))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zmc[,i] <- Zmc[,i-1]*((1-tau)*zcm)/(i-1)
    }
  }

  ## Compute Laguerre Polynomials
  Lpc <- Zpc%*%t(B)
  Lmc <- Zmc%*%t(B_tilde)

  ## Compute the density
  fc <- rep(0,length(z_censored))
  fc[z_censored>=0] <- ((1-tau)/(scp))*tau*exp(-tau*zcp)*(Lpc%*%theta)^2
  fc[z_censored< 0] <- (tau*(1-tau)/scn)*exp((1-tau)*zcm)*(Lmc%*%theta_tilde)^2

  ## Derivative of log-distribution
  logFgrad <- Her[Delta==0,]*((fc/(1-F))*(z[Delta==0]/sigma[Delta==0])*dsigma[Delta==0])

  ## Derivative of log-likelihood
  Lgrad_sigma <- colSums(logfgradp)+colSums(logfgradm)+colSums(logFgrad)

  if(derivative=="sigma") {
    ## Return
    return(list("objective"=L/n,"gradient"=Lgrad_sigma/n))
  }


  if(derivative=="theta" | derivative=="all" | derivative=="beta_theta") {
    #### Compute the derivative of L with respect to theta
    ##  Derivative of f with respect to theta
    fgrad <- 2*(1/sup)*(1-tau)*tau*as.numeric(exp(-tau*zup)*(Lp%*%theta))*Lp
    ## of F with respect to theta
    Fgrad <- rep(0,m+1)
    for(i1 in 0:m) {
      for(i2 in 0:m) {
        Fgrad <- Fgrad+(1-tau)*v[i1+1]*v[i2+1]*as.matrix(Ilower[,i1+i2+1])%*%t(sum(theta*B[,i2+1])*B[,i1+1]+sum(theta*B[,i1+1])*B[,i2+1])
      }
    }

    ## Compute gradient of likelihood
    Lgrad_theta <- colSums(fgrad/f[z_uncensored>=0])-colSums(Fgrad/(1-F[z_censored>=0]))

    if(derivative=="theta") {
      ## Return
      return(list("objective"=L/n,"gradient"=Lgrad_theta/n))
    }
  }
  if(derivative=="theta_tilde" | derivative=="all" | derivative=="beta_theta_tilde") {
    #### Compute the derivative of L with respect to theta_tilde
    ##  Derivative of f with respect to theta_tilde
    fgrad <- 2*(1/sun)*(1-tau)*tau*as.numeric(exp((1-tau)*zum)*(Lm%*%theta_tilde))*Lm
    ## of F with respect to theta_tilde
    Fgrad <- rep(0,m_tilde+1)
    for(i1 in 0:m_tilde) {
      for(i2 in 0:m_tilde) {
        Fgrad <- Fgrad+tau*v_tilde[i1+1]*v_tilde[i2+1]*as.matrix(Iupper[,i1+i2+1])%*%t(sum(theta_tilde*B_tilde[,i2+1])*B_tilde[,i1+1]+sum(theta_tilde*B_tilde[,i1+1])*B_tilde[,i2+1])
      }
    }

    ## Compute gradient of likelihood
    Lgrad_theta_tilde <- colSums(fgrad/f[z_uncensored<0])-colSums(Fgrad/(1-F[z_censored<0]))

    if(derivative=="theta_tilde") {
      ## Return
      return(list("objective"=L/n,"gradient"=Lgrad_theta_tilde/n))
    }
  }
  if(derivative=="beta" | derivative=="all" | derivative=="beta_theta" | derivative=="beta_theta_tilde") {
    ## Derivative of log-density with respect to beta
    if(sum(up_index)>=1) {
      if(m==0) {
        logfgradp <- (X[up_index,])*(tau/sup)
      } else {
        if(m==1) {
          helper <-  t(t(Zp[,1:m      ]))
        } else if(dim(Zp)[1]==1) {
          helper <-  t(Zp[,1:m      ])
        } else {
          helper <-  Zp[,1:m      ]
        }
        logfgradp <-     (X[up_index,])*(tau*as.numeric(1+2*helper%*%(t(B      )%*%theta      )[2:(m      +1)]/(Lp%*%theta      ))/sup)
      }
      if(sum(up_index)==1) {
        logfgradp <- matrix(logfgradp,nrow=1)
      }
    } else {
      logfgradp <- matrix(0,nrow=1,ncol=length(beta))
    }

    if(sum(um_index)>=1) {
      if(m_tilde==0) {
        logfgradm <- (X[um_index,])*(-(1-tau)/sun)
      } else {
        if(m_tilde==1) {
          helper <- t(t(Zm[,1:m_tilde]))
        } else if(dim(Zm)[1]==1) {
          helper <- t(Zm[,1:m_tilde])
        } else {
          helper <- Zm[,1:m_tilde]
        }
        logfgradm <- (X[um_index,])*(-(1-tau)*as.numeric(1+2*helper%*%(t(B_tilde)%*%theta_tilde)[2:(m_tilde+1)]/(Lm%*%theta_tilde))/sun)
      }
      if(sum(um_index)==1) {
        logfgradm <- matrix(logfgradm,nrow=1)
      }
    } else {
      logfgradm <- matrix(0,nrow=1,ncol=length(beta))
    }

    ## Compute the density for the censored observations
    ## Compute the powers of z
    Zpc <- matrix(1,ncol=m+1,nrow=length(zcp))
    if(m+1>=2) {
      for(i in 2:(m+1)) {
        Zpc[,i] <- Zpc[,i-1]*(-tau*zcp)/(i-1)
      }
    }
    Zmc <- matrix(1,ncol=m_tilde+1,nrow=length(zcm))
    if(m_tilde+1>=2) {
      for(i in 2:(m_tilde+1)) {
        Zmc[,i] <- Zmc[,i-1]*((1-tau)*zcm)/(i-1)
      }
    }

    ## Compute Laguerre Polynomials
    Lpc <- Zpc%*%t(B)
    Lmc <- Zmc%*%t(B_tilde)

    ## Compute the density
    fc <- rep(0,length(z_censored))
    fc[z_censored>=0] <- ((1-tau)/(scp))*tau*exp(-tau*zcp)*(Lpc%*%theta)^2
    fc[z_censored< 0] <- (tau*(1-tau)/scn)*exp((1-tau)*zcm)*(Lmc%*%theta_tilde)^2

    ## Derivative of log-distribution
    logFgrad <- X[Delta==0,] * fc/(1-F)

    ## Derivative of log-likelihood
    Lgrad_beta <- colSums(logfgradp)+colSums(logfgradm)+colSums(logFgrad)

    if(derivative=="beta") {
      ## Return
      return(list("objective"=L/n,"gradient"=c(Lgrad_beta,Lgrad_sigma)/n))
    }
  }

  if(derivative=="beta_theta") {
    return(list("objective"=L/n,"gradient"=c(Lgrad_beta,Lgrad_theta,Lgrad_sigma)/n))
  }
  if(derivative=="beta_theta_tilde") {
    return(list("objective"=L/n,"gradient"=c(Lgrad_beta,Lgrad_theta_tilde,Lgrad_sigma)/n))
  }
  return(list("objective"=L/n,"gradient"=c(Lgrad_beta,Lgrad_theta,Lgrad_theta_tilde,Lgrad_sigma)/n))
}


## Computes the log-likelihood of the data for a given set of parameters in the
## same way as likelihood above does but theta and theta_tilde are provided in
## polar coordinates.
## Input:
##  tau         - see likelihood
##  beta        - see likelihood
##  X           - see likelihood
##  Y           - see likelihood
##  Delta       - see likelihood
##  derivative  - see likelihood
##  theta_polar - Contains the angles of the polar coordinate representation of
##                theta, the radius must always be equal to one. If theta=1
##                shall be provided, set here theta_tilde=FALSE.
##  theta_tilde_polar - Same for theta_tilde as theta_polar for theta
## Output: See likelihood
likelihood_polar2 <- function(tau,beta,theta_polar,theta_tilde_polar,H_polar,X,X_s,type,Y,Delta,derivative=FALSE, link="exp",deg) {
  if(isFALSE(theta_polar)==TRUE) {
    theta <- 1
    if(derivative=="all") {
      derivative <- "beta_theta_tilde"
    }
  } else {
    theta <- SphericalCubature::polar2rect(1,theta_polar)
  }
  if(isFALSE(theta_tilde_polar)==TRUE) {
    theta_tilde <- 1
    if(derivative=="all") {
      derivative <- "beta_theta"
    }
  } else {
    theta_tilde <- SphericalCubature::polar2rect(1,theta_tilde_polar)
  }

  H = SphericalCubature::polar2rect(1, H_polar)

  out <- likelihood2(tau,beta,theta,theta_tilde,H,X,X_s,type,Y,Delta,derivative=derivative, link=link,deg=deg)

  if(derivative==FALSE) {
    return(out)
  }
  if(derivative=="theta") {
    DT <- polar_derivative(theta_polar)
    DT_H = polar_derivative(H_polar)
    grad <- t(DT)%*%out$gradient
    grad_H = t(DT_H)%*%out$gradient
    return(list("objective"=out$objective,"gradient"=c(grad, grad_H)))
  }
  if(derivative=="theta_tilde") {
    DT <- polar_derivative(theta_tilde_polar)
    DT_H = polar_derivative(H_polar)
    grad <- t(DT)%*%out$gradient
    grad_H = t(DT_H)%*%out$gradient
    return(list("objective"=out$objective,"gradient"=c(grad, grad_H)))
  }
  if(derivative=="beta") {
    DT_H = polar_derivative(H_polar)
    p1 = length(beta)
    p2 = length(H)
    grad = rep(0, p1 + p2 - 1)
    grad[1:p1] = out$gradient[1:p1]
    grad[(p1+1):(p1+p2-1)] = t(DT_H)%*%out$gradient[(p1+1):(p1+p2)]
    return(list("objective"=out$objective, "gradient"=grad))
  }
  if(derivative=="all") {
    DT <- polar_derivative(theta_polar)
    DT_tilde <- polar_derivative(theta_tilde_polar)
    DT_H = polar_derivative(H_polar)
    p1 <- length(beta)
    p2 <- length(theta)
    p3 <- length(theta_tilde)
    p4 = length(H)
    grad <- rep(0,p1+p2+p3+p4-3)
    grad[1:p1] <- out$gradient[1:p1]
    grad[(p1+1 ):(p1+p2   -1)] <- t(DT      )%*%out$gradient[(p1   +1):(p1+p2)]
    grad[(p1+p2):(p1+p2+p3-2)] <- t(DT_tilde)%*%out$gradient[(p1+p2+1):(p1+p2+p3)]
    grad[(p1+p2+p3-1):(p1+p2+p3+p4-3)] <- t(DT_H)%*%out$gradient[(p1+p2+p3+1):(p1+p2+p3+p4)]
    return(list("objective"=out$objective,"gradient"=grad))
  }
  if(derivative=="beta_theta") {
    DT <- polar_derivative(theta_polar)
    DT_H = polar_derivative(H_polar)
    p1 <- length(beta)
    p2 <- length(theta)
    p3= length(H)
    grad <- rep(0,p1+p2+p3-2)
    grad[1:p1] <- out$gradient[1:p1]
    grad[(p1+1 ):(p1+p2   -1)] <- t(DT)%*%out$gradient[(p1+1):(p1+p2)]
    grad[(p1+p2 ):(p1+p2  +p3 -2)] <- t(DT_H)%*%out$gradient[(p1+p2+1):(p1+p2+p3)]
    return(list("objective"=out$objective,"gradient"=grad))
  }
  if(derivative=="beta_theta_tilde") {
    DT_tilde <- polar_derivative(theta_tilde_polar)
    DT_H = polar_derivative(H_polar)
    p1 <- length(beta)
    p2 <- length(theta_tilde)
    p3= length(H)
    grad <- rep(0,p1+p2+p3-2)
    grad[1:p1] <- out$gradient[1:p1]
    grad[(p1+1):(p1+p2-1)] <- t(DT_tilde)%*%out$gradient[(p1+1):(p1+p2)]
    grad[(p1+p2 ):(p1+p2  +p3 -2)] <- t(DT_H)%*%out$gradient[(p1+p2+1):(p1+p2+p3)]
    return(list("objective"=out$objective,"gradient"=grad))
  }
}
