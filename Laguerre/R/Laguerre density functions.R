#' Computes f_{theta,theta_tilde}(z), i.e.,  the Laguerre density as specified
#' in the paper.
#' Input:
#'  z           - Vector of points at which the density shall be evaluated
#'  tau,theta,  - see likelihood
#'  theta_tilde
#' Output: A vector of the same length as z which contains the corresponding
#'         values of the density with quantile tau and parameters theta and
#'         theta_tilde
#' @export
laguerre_density <- function(z,tau,theta,theta_tilde) {
  ## Prepare Data matrices
  m <- length(theta)-1
  m_tilde <- length(theta_tilde)-1
  n <- length(z)

  zp <- z[z>=0]
  zm <- z[z <0]

  ## Compute the matrix of binomial coefficients
  B       <- matrix(0,ncol=m+1      ,nrow=m+1)
  B_tilde <- matrix(0,ncol=m_tilde+1,nrow=m_tilde+1)

  for(k in 1:(m+1)) {
    B[k,1:k] <- choose(k-1,0:(k-1))
  }
  for(k in 1:(m_tilde+1)) {
    B_tilde[k,1:k] <- choose(k-1,0:(k-1))
  }

  ## Compute the powers of z
  Zp <- matrix(1,ncol=m+1,nrow=length(zp))
  if(m+1>=2) {
    for(i in 2:(m+1)) {
      Zp[,i] <- Zp[,i-1]*(-tau*zp)/(i-1)
    }
  }
  Zm <- matrix(1,ncol=m_tilde+1,nrow=length(zm))
  if(m_tilde+1>=2) {
    for(i in 2:(m_tilde+1)) {
      Zm[,i] <- Zm[,i-1]*((1-tau)*zm)/(i-1)
    }
  }

  ## Compute Laguerre Polynomials
  Lp <- Zp%*%t(B)
  Lm <- Zm%*%t(B_tilde)

  ## Compute the density
  f <- rep(0,n)
  f[z>=0] <- (1-tau)*tau*exp(-tau*zp)*(Lp%*%theta)^2
  f[z< 0] <- tau*(1-tau)*exp((1-tau)*zm)*(Lm%*%theta_tilde)^2

  return(f)
}




###### FUNCTIONS TO CALCULATE THE N-TH MOMENT OF THE LAGUERRE DENSITY #########

### Gamma term

gamm = function(n, i, k , r, s){
  gamm = (gamma(n+r+s+1)*choose(i, r)*choose(k, s)*(-1)^(r+s))/(factorial(r)*factorial(s))
  return(gamm)
}


### Function for I_ik integral
I_ik = function(n,i,k){
  terms = matrix(0, nrow= i+1, ncol = k+1)
  for (a in 0:i){
    for(b in 0:k){
      terms[a+1,b+1] = gamm(n, i, k, a, b)
    }
  }

  result = sum(colSums(terms))
  return(list(result, terms))
}

#Moments of the Laguerre density
#' @export
laguerre_moment = function(n, theta, theta_tilde, tau){
  terms_theta = matrix(0, ncol=length(theta), nrow =length(theta) )
  terms_theta_tilde = matrix(0, ncol = length(theta_tilde), nrow = length(theta_tilde))

  for (a in 1:length(theta_tilde)){
    for (b in 1:length(theta_tilde)){
      terms_theta_tilde[a, b] = theta_tilde[a]*theta_tilde[b]*I_ik(n,a-1,b-1)[[1]]
    }
  }

  for (a in 1:length(theta)){
    for (b in 1:length(theta)){
      terms_theta[a, b] = theta[a]*theta[b]*I_ik(n,a-1,b-1)[[1]]
    }
  }

  result = (((-1)^n)*tau/(1-tau)^n)*sum(colSums(terms_theta_tilde)) +
    ((1-tau)/tau^n)*sum(colSums(terms_theta))

  return(result)

}

#Variance of the Laguerre density
#' @export
laguerre_var = function(theta, theta_tilde, tau){
  var = laguerre_moment(2,theta,theta_tilde,tau) - laguerre_moment(1,theta,theta_tilde,tau)^2
  return(var)
}
