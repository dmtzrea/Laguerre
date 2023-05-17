## Wrapper functions which make the syntax more suitable for the optimization
## and return the negative log-likelihood.
likelihood_wrapper_beta_only2 <- function(x,X,X_s,type,Y,Delta,tau, dims, link="exp",deg) {
  beta <- x[1:dims[1]]
  H_polar = x[(dims[1]+1):(dims[1]+dims[4])]
  out <- likelihood_polar2(tau,beta,FALSE,FALSE,H_polar,X,X_s,type,Y,Delta,derivative="beta", link=link,deg=deg)
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}

likelihood_wrapper_beta_theta2 <- function(x,X,X_s,type,Y,Delta,tau,dims,link="exp",deg) {
  beta <- x[1:dims[1]]
  theta_polar <- x[(dims[1]+1):(dims[1]+dims[2])]
  H_polar = x[(dims[1]+dims[2]+dims[3]+1):(dims[1]+dims[2]+dims[3]+dims[4])]
  out <- likelihood_polar2(tau,beta,theta_polar,FALSE,H_polar,X,X_s,type,Y,Delta,derivative="all", link=link,deg=deg)
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}

likelihood_wrapper_beta_theta_tilde2 <- function(x,X,X_s,type,Y,Delta,tau,dims, link="exp",deg) {
  beta <- x[1:dims[1]]
  theta_tilde_polar <- x[(dims[1]+dims[2]+1):(dims[1]+dims[2]+dims[3])]
  H_polar = x[(dims[1]+dims[2]+dims[3]+1):(dims[1]+dims[2]+dims[3]+dims[4])]
  out <- likelihood_polar2(tau,beta,FALSE,theta_tilde_polar,H_polar,X,X_s,type,Y,Delta,derivative="beta_theta_tilde", link=link,deg=deg)
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}

likelihood_wrapper_all2 <- function(x,X,X_s,type,Y,Delta,tau,dims,link="exp",deg) {
  beta <- x[1:dims[1]]
  theta_polar <- x[(dims[1]+1):(dims[1]+dims[2])]
  theta_tilde_polar <- x[(dims[1]+dims[2]+1):(dims[1]+dims[2]+dims[3])]
  H_polar = x[(dims[1]+dims[2]+dims[3]+1):(dims[1]+dims[2]+dims[3]+dims[4])]
  out <- likelihood_polar2(tau,beta,theta_polar,theta_tilde_polar,H_polar, X,X_s,type,Y,Delta,derivative="all", link=link,deg=deg)
  return(list("objective"=-out$objective,"gradient"=-out$gradient))
}
