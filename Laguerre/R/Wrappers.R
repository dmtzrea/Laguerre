## Wrapper functions which make the syntax more suitable for the optimization
## and return the negative log-likelihood.
likelihood_wrapper_beta_only <- function(beta,X,Y,Delta,tau) {
  out <- likelihood(tau,beta,1,1,X,Y,Delta,derivative="beta")
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}
likelihood_wrapper_beta_theta <- function(par,X,Y,Delta,tau,m) {
  p <- length(par)-m
  out <- likelihood_polar(tau,par[1:p],par[(p+1):(p+m)],FALSE,X,Y,Delta,derivative="all")
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}
likelihood_wrapper_beta_theta_tilde <- function(par,X,Y,Delta,tau,m_tilde) {
  p <- length(par)-m_tilde
  out <- likelihood_polar(tau,par[1:p],FALSE,par[(p+1):(p+m_tilde)],X,Y,Delta,derivative="all")
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}
likelihood_wrapper_all <- function(x,X,Y,Delta,tau,dims) {
  beta <- x[1:dims[1]]
  theta_polar <- x[(dims[1]+1):(dims[1]+dims[2])]
  theta_tilde_polar <- x[(dims[1]+dims[2]+1):(dims[1]+dims[2]+dims[3])]

  out <- likelihood_polar(tau,beta,theta_polar,theta_tilde_polar,X,Y,Delta,derivative="all")
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}
