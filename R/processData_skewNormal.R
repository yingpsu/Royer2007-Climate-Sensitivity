##==============================================================================
## processData_skewNormal.R
##
## Read CO2 proxy data, with age/amount uncertainties.
##
## Fit skew-normal distributions at each time point, assuming the +sigma1/-sigma2
## range contains 0.68 probability mass. Want to fit the 16%, 84% quantiles and
## the mean. Skew-normal mean is xi+omega*delta*sqrt(2/pi), where
## delta=alpha/sqrt(1+alpha^2).
##
## Write out a file with the skew-normal distribution parameters on it, to use
## in the calibration.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

rm(list=ls())

library(sn)
library(DEoptim)

# col 1 is age (million years ago), col 2 is CO2 (ppmv), col 3 is lower error
# bound, col 4 is upper error bound
dat <- read.csv('../input_data/CO2_Proxy_Foster2017_calib.csv', fill=TRUE, header=TRUE)

# want to determine the three skew-normal parameters (xi, omega, alpha) that
# give the median and 68% probability mass between upper/lower erorr bounds
# (dat$co2_low and dat$co2_high, with dat$co2 as the mean)

# define a function for fitting the parameters within a given tolerance
# parameters[1]=xi (location)
# parameters[2]=omega (scale)
# parameters[3]=alpha (slant)
sse_quantiles <- function(parameters, q16, q50, q84) {
    q16.fit <- qsn(0.16, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    q84.fit <- qsn(0.84, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    # fit based on median
    #q50.fit <- qsn(0.50, xi=parameters[1], omega=parameters[2], alpha=parameters[3], tol=1e-8, solver="NR")
    # fit based on mean
    q50.fit <- parameters[1] + parameters[2]*parameters[3]*sqrt(2/(pi*1+parameters[3]^2))
    sse <- sqrt((q16.fit-q16)^2 + (q50.fit-q50)^2 + (q84.fit-q84)^2)
    return(sse)
}

# preliminary testing to see how many iterations are needed to get convergence
bound.lower <- c(-1e4, 0, -100)
bound.upper <- c(1e4, 1e4, 100)

# testing
niter.test <- c(100,200,300,400,500,600,700,800,900,1000)#,2000,5000)
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
                      q16=dat$co2_low[1], q50=dat$co2[1], q84=dat$co2_high[1])
  t2 <- proc.time()
  parameters <- outDEoptim$optim$bestmem
  error.test[i] <- sse_quantiles(parameters, dat$co2_low[1], dat$co2[1], dat$co2_high[1])
  times.test[i] <- t2[3]-t1[3]
}

cbind(niter.test, error.test, times.test)

# for the first data point, >=300 iterations stably gives ~0 error in the
# fitted quantiles. this should be more than good enough.

# fit the CO2 data
niter.deoptim <- 400
parameters.co2 <- vector('list',3)
names(parameters.co2) <- c('xi','omega','alpha')
for (p in 1:3) {parameters.co2[[p]] <- rep(NA,length(dat$co2))}
error.co2 <- rep(NA,length(dat$co2))
t1 <- proc.time()
pb <- txtProgressBar(min=0,max=length(dat$co2),initial=0,style=3)
for (i in 1:length(dat$co2)){
  outDEoptim <- DEoptim(sse_quantiles, bound.lower, bound.upper,
                      DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,
                      CR=CR.deoptim,trace=FALSE),
                      q16=dat$co2_low[i], q50=dat$co2[i], q84=dat$co2_high[i])
  parameters.co2$xi[i] <- outDEoptim$optim$bestmem[1]
  parameters.co2$omega[i] <- outDEoptim$optim$bestmem[2]
  parameters.co2$alpha[i] <- outDEoptim$optim$bestmem[3]
  error.co2[i] <- sse_quantiles(outDEoptim$optim$bestmem, dat$co2_low[i], dat$co2[i], dat$co2_high[i])
  setTxtProgressBar(pb, i)
}
close(pb)
t2 <- proc.time()

# save workspace image
save.image(file = "fit_co2_data_skewNormal.RData")

# write new CSV file, copy of Foster calib data set, but with the three
# skew-normal parameters for each CO2 data point.
dat.co2 <- cbind(dat,parameters.co2$xi, parameters.co2$omega, parameters.co2$alpha)
colnames(dat.co2) <- c(colnames(dat),'xi_co2','omega_co2','alpha_co2')

# clean out the rows that have co2 <= 0, or co2_low=co2_high=0. there are not useful.
ind_co2 <- which(dat.co2$co2 <= 0)   # for default experiments
ind_co2_low <- which(dat.co2$co2_low==0)
ind_co2_high <- which(dat.co2$co2_high==0)
irem <- unique(c( intersect(ind_co2_low,ind_co2_high),ind_co2))
dat.co2 <- dat.co2[-irem,]

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.out <- paste('../input_data/CO2_Proxy_Foster2017_calib_SN-co2_',today,'.csv',sep='')
write.csv(dat.co2, file=filename.out)

##==============================================================================
## End
##==============================================================================
