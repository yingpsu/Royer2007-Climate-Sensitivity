##==============================================================================
## Read CO2 proxy data
## Fit skew-normal distributions at each time point, assuming the +sigma1/-sigma2
## range contains 0.68 probability mass.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

library(sn)

# col 1 is age (million years ago), col 2 is CO2 (units?), col 3 is lower error
# bound, col 4 is upper error bound
dat <- read.table('../data/proxydata_clean.txt', fill=TRUE, header=TRUE)

# want to determine the three skew-normal parameters (xi, omega, alpha) that
# give the median and 68% probability mass between upper/lower erorr bounds.

# define a function for fitting the parameters within a given tolerance
# parameters[1]=xi (location)
# parameters[2]=omega (scale)
# parameters[3]=alpha (slant)
sse_quantiles <- function(parameters, q16, q50, q84) {
    q16.fit <- qsn(0.16, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    q50.fit <- qsn(0.50, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    q84.fit <- qsn(0.84, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    sse <- sqrt((q16.fit-q16)^2 + (q50.fit-q50)^2 + (q84.fit-q84)^2)
    return(sse)
}

# preliminary testing to see how many iterations are needed to get convergence
bound.lower <- c(0, 0, -1e4)
bound.upper <- c(1e4, 1e4, 100)

# testing
niter.test <- c(100,200,500,1000,2000,5000)
error.test <- rep(NA, length(niter.test))
times.test <- rep(NA, length(niter.test))
#niter.deoptim=500        # number of iterations for DE optimization
NP.deoptim=30             # population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8             # as suggested by Storn et al (2006)
CR.deoptim=0.9            # as suggested by Storn et al (2006)

for (i in 1:length(niter.test)){
  t1 <- proc.time()
  outDEoptim <- DEoptim(sse_quantiles, bound.lower, bound.upper,
                      DEoptim.control(NP=NP.deoptim,itermax=niter.test[i],F=F.deoptim,
                      CR=CR.deoptim,trace=FALSE),
                      q16=dat$unc_low[1], q50=dat$CO2[1], q84=dat$unc_high[1])
  t2 <- proc.time()
  parameters <- outDEoptim$optim$bestmem
  error.test[i] <- sse_quantiles(parameters, dat$unc_low[1], dat$CO2[1], dat$unc_high[1])
  times.test[i] <- t2[3]-t1[3]
}


# TODO:
# 1. multivariate (both CO2 and age uncertainties) (qsn(...) -> qmsn(...))
# 2. fit for each point in time

##==============================================================================
## End
##==============================================================================
