#' A function that plots the estimated quantile and sigma.
#' The arguments are:
#' H: The estimated coefficients of the Hermite series for sigma.
#' beta: The estimated regression coefficients.
#' X: The design matrix of the regression.
#' Y: The vector of responses.
#' Delta: The censorship indicator.
#' @export
graph = function(fit){
  H = fit$H
  beta = fit$beta
  X = fit$X
  X_s = fit$X_s
  Y = fit$Y
  link = fit$link
  type = fit$type
  deg = fit$deg
  Her = Her(X_s, deg=deg, type=type)
  if(link == "exp"){
    sigma = exp(Her%*%H)/exp(1)
  }

  if(link=="quad"){
    sigma = (Her%*%H)^2
  }

  if (link!="exp" & link!="quad"){
    sigma = as.vector(unlist(lapply(link, function(f) f(Her%*%H))[1]))
    dsigma = as.vector(unlist(lapply(link, function(f) f(Her%*%H))[2]))
  }
  sigma = sigma*sqrt(laguerre_var(fit$theta, fit$theta_tilde,fit$quantile))
  s = as.matrix(seq(from=min(X[,2]), to=max(X[,2]), by=0.01))
  s_smooth = (exp(Her(s, deg = deg, type=type)%*%H)/exp(1))*sqrt(laguerre_var(fit$theta, fit$theta_tilde,fit$quantile))
  data = as.data.frame(cbind(X[,2],X_s, sigma,X%*%beta, Y,Delta))
  colnames(data) = c("x", "x_s", "sigma", "quantile", "Y", "Delta")
  data2 = as.data.frame(cbind(s,s_smooth))
  colnames(data2) = c("s", "s_smooth")
  if(dim(X_s)[2]==1){
    print(ggplot2::ggplot(data, ggplot2::aes_string(x = "x", y = "sigma")) + ggplot2::geom_point() + ggplot2::geom_line(color="orange") +
            ggplot2::ggtitle("sigma(x)"))

    print(ggplot2::ggplot(data2, ggplot2::aes_string(x = "s", y = "s_smooth")) + ggplot2::geom_point() + ggplot2::geom_line(color="green") +
            ggplot2::ggtitle("sigma(x) smooth"))
  }
  print(ggplot2::ggplot(data, ggplot2::aes_string(x = "x", y = "Y", colour = as.factor(Delta))) + ggplot2::geom_point() +
          ggplot2::geom_line(data=data, ggplot2::aes_string(x = "x", y = "quantile"), color="black") +
          ggplot2::ggtitle("Regression Quantile"))
}
