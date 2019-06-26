##==============================================================================
## dasn.R
##
## Approximate skew-normal, as in Ashour and Abdel-hameed [2010]:
## https://www.researchgate.net/publication/232415595_Approximate_skew_normal_distribution
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


##==============================================================================
## testing...
##==============================================================================

## convolution?

lb <- -50; ub <- 200
x = seq(lb,ub,1)
m1 = 40; y1 = dnorm(x, mean=m1, sd=20); s1 = rnorm(n=10000, mean=m1, sd=20)
m2 = 0; y2 = dnorm(x, mean=m2, sd=10); s2 = rnorm(n=10000, mean=m2, sd=10); ss = s1+s2; y3 = density(ss, from=lb, to=ub)
plot(y3, type='l', ylim=c(0,.08), main=''); lines(x,y1, type='l', lty=2); lines(x,y2,type='l', lty=2);
z = convolve(y1,y2); lines(z, col='purple')





co2 <- seq(from=-10000,to=15000,by=10)

# a representative cross-sectional plot
t <- 34
f_fit <- likelihood_fit[[t]](co2)
f_fit[is.na(f_fit)] <- 0
plot(co2,f_fit*10, xlab='CO2 (ppmv)', ylab='density', type='l', ylim=c(0,0.006))

# add in normal errors
unc <- dnorm(co2, mean=0, sd=100); lines(co2,unc*10, col='purple')
z <- convolve(f_fit*10, unc*10); lines(co2, z, lty=2)




# looks good:

dx <- 10
co2 <- seq(from=-10000,to=15000,by=dx)
t <- 35
sigma <- 400

# ONE WAY: (do integration within `convolve`)

f_fit <- likelihood_fit[[t]](co2)
f_fit[is.na(f_fit)] <- 0

unc <- dnorm(co2, mean=0, sd=sigma);
z <- convolve(f_fit*dx, unc*dx);
i0 <- which(co2==0); z <- c(z[(length(z)-i0):length(z)], z[1:(length(z)-i0-1)])

ll[t] <- approx(co2, z, xout=model_out[,"co2"], method = "linear", yleft=0, yright=0, rule=2, ties = mean)


# ANOTHER WAY: (explicitly do the integration)
#f.like <- function(co2) {idx <- which(co2 >= 0 & co2 <= 10000); output <- rep(0, length(co2)); output[idx] <- likelihood_fit[[t]](x[idx]); return(output)};        # normal (mu=1.5, sigma=0.5)
f.Y <- function(y) dnorm(y,0,sigma)   # log-normal (mu=1.5, sigma=0.75)
# convolution integral
f.Z <- function(z) integrate(function(x,z) f.Y(z-x)*f.X(x),-2000,10000,z)$value
f.Z <- Vectorize(f.Z)                    # need to vectorize the resulting fn.


plot(co2,f.Z(co2),type='l',lty=2,col="red", ylim=c(0,0.0007)); lines(co2,f.X(co2))

#plot(co2,f_fit*dx, xlab='CO2 (ppmv)', ylab='density', type='l', ylim=c(0,0.006))
#lines(co2,unc*dx, col='purple')
#lines(co2, z, lty=2)


##==============================================================================
## end testing
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
