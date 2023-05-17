### A function that calculates the cross products of orthogonal polynomials
### up to a certain order.
#' @importFrom polynom polynomial
#' @export
Her = function(X_s, deg, type){
  ### Computing Hermite polynomials
  ### These are the terms in each variable up to the specified degree, without cross-products
  Her = c()
  for(i in 1:dim(X_s)[2]){
    Her1 = polys(X_s[,i], which = type[i], normalized = TRUE, n = deg)
    #for (k in 1:(deg+1)){
    #Remember the factor of 5 here
    #Her1[, k] = 5*EQL::hermite(X_s[,i], k-1, prob=FALSE)/sqrt((sqrt(pi)*2^(k-1)*factorial(k-1)))
    #}

    Her = cbind(Her, Her1)
  }

  ###Generating cross-products
  if(dim(X_s)[2]>1){
    Her2 = c()
    for(j in 1:dim(X_s)[1]){
      prov = Her[j, 1:(1+deg)]
      for (i in seq(from=1+deg+1,to=(dim(Her)[2]-deg), by = deg+1)){
        prov = prov%o%Her[j, i:(i+deg)]
      }
      ind = as.matrix((which(prov == prov,arr.ind = TRUE) -1))
      ind = ind[which(rowSums(ind)<=deg),]
      prov = prov[as.matrix(ind)+1]
      Her2=rbind(Her2,prov)
    }
    Her = Her2
  }


  return(Her)
}


### Function to evaluate the orthogonal polynomial
###

polys = function(x, which, normalized = TRUE, n){
  if(which == "Laguerre"){
    leg4coef <- orthopolynom::laguerre.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }
  if(which == "Hermite"){
    leg4coef <- orthopolynom::hermite.h.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }
  if(which == "Legendre"){
    leg4coef <- orthopolynom::legendre.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=scaleX(x,u=-1,v=1))))
  }

  if(which == "orth"){
    leg4 <- cbind(rep(1,length(x)), stats::poly(x,degree=n))
  }

  colnames(leg4) = 1:(n+1)
  return(leg4)
}

