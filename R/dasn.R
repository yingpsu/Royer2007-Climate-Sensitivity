##==============================================================================
## dasn.R
##
## Approximate skew-normal, as in Ashour and Abdel-hameed [2010]:
## https://www.researchgate.net/publication/232415595_Approximate_skew_normal_distribution
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

##==============================================================================
dasn_single <- function(z, phi, alpha){

    if (alpha >= 0) {
        if (z < -3/alpha) {
            return(0)
        } else if (z < -1/alpha) {
            return((phi/8)*(9*alpha*z + 3*(alpha*z)^2 + ((alpha*z)^3)/3 + 9))
        } else if (z < 1/alpha) {
            return((phi/4)*(3*alpha*z - ((alpha*z)^3)/3 + 4))
        } else if (z < 3/alpha) {
            return((phi/8)*(9*alpha*z - 3*(alpha*z)^2 + ((alpha*z)^3)/3 + 7))
        } else {
            return(2*phi)
        }
    } else {
        if (z > -3/alpha) {
            return(0)
        } else if (z > -1/alpha) {
            return((phi/8)*(9*alpha*z + 3*(alpha*z)^2 + ((alpha*z)^3)/3 + 9))
        } else if (z > 1/alpha) {
            return((phi/4)*(3*alpha*z - ((alpha*z)^3)/3 + 4))
        } else if (z > 3/alpha) {
            return((phi/8)*(9*alpha*z - 3*(alpha*z)^2 + ((alpha*z)^3)/3 + 7))
        } else {
            return(2*phi)
        }
    }
}
##==============================================================================


dasn <- function(x, xi, omega, alpha){

  z <- (x-xi)/omega                         # transform to center and scale
  phi <- exp(-z*z/2)/sqrt(2*pi*omega*omega) # constant needed in all cases

  return(sapply(1:length(z), FUN=function(i){dasn_single(z[i], phi[i], alpha)}))
}

##==============================================================================
if(FALSE) {
dasn <- function(x, xi, omega, alpha){

    z <- (x-xi)/omega             # transform to center and scale
    phi <- exp(-z*z/2)/sqrt(2*pi*omega*omega) # constant needed in all cases
    f_out <- rep(0, length(z))    # initialize output

    # if (z < -3/alpha)
    idx <- which(z < -3/alpha)
    if(length(idx) > 0) f_out[idx] <- 0

    # if (z < -1/alpha)
    idx <- which(z < -1/alpha & z >= -3/alpha)
    if(length(idx) > 0) f_out[idx] <- (phi[idx]/8)*(9*alpha*z[idx] + 3*(alpha*z[idx])^2 + ((alpha*z[idx])^3)/3 + 9)

    # if (z < 1/alpha)
    idx <- which(z < 1/alpha & z >= -1/alpha)
    if(length(idx) > 0) f_out[idx] <- (phi[idx]/4)*(3*alpha*z[idx] - ((alpha*z[idx])^3)/3 + 4)

    # if (z < 3/alpha)
    idx <- which(z < 3/alpha & z >= 1/alpha)
    if(length(idx) > 0) f_out[idx] <- (phi[idx]/8)*(9*alpha*z[idx] - 3*(alpha*z[idx])^2 + ((alpha*z[idx])^3)/3 + 7)

    # if (z >= 3/alpha)
    idx <- which(z >= 3/alpha)
    if(length(idx) > 0) f_out[idx] <- 2*phi[idx]

    return(f_out)
}
}
