################################################################################
## The functions below are used internally in the functions above. In most    ##
## cases the user will not need to call these functions directly.             ##
################################################################################

## Computes the derivative of the polar coordinate transform for r=1 and set of
## angles phi
polar_derivative <- function(phi) {
  m <- length(phi)
  deriv <- matrix(0,ncol=m,nrow=m+1)
  for(i in 1:m) {
    shift <- rep(0,m)
    shift[i] <- pi/2
    deriv[i:(m+1),i] <- SphericalCubature::polar2rect(1,phi+shift)[i:(m+1)]
  }

  return(deriv)
}
