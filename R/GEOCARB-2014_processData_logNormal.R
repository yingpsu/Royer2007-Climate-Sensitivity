##==============================================================================
## GEOCARB-2014_processData_logNormal.R
##
## Read CO2 proxy data, with age/amount uncertainties.
##
## Fit log-normal distributions at each time point, assuming the +sigma1/-sigma2
## range contains 0.68 probability mass. The fit is done using evolutionary
## algorithm for optimization. The default general-purpose `optim` routine in R
## yielded poor results, so `DEoptim` is used instead
##
## Write out a file with the gamma distribution parameters on it, to use
## in the calibration.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

library(DEoptim)

# col 1 is age (million years ago), col 2 is CO2 (units?), col 3 is lower error
# bound, col 4 is upper error bound
dat <- read.csv('../input_data/CO2_Proxy_Foster2017_calib.csv', fill=TRUE, header=TRUE)

# first, eliminate all of the pCO2 <= 0 points
ind_remove <- which(dat$co2 <= 0)
dat <- dat[-ind_remove,]

# now, deal with the data point with co2_low = 0. we will fit the 16% quantile
# to the co2_low value, so it needs to be away from 0. Also, 0 is not a proper
# lower uncertainty bound unless infinity is the appropriate upper bound. So,
# use the log-normal rule:
# co2_high = co2*s --> s = co2_high/co2 --> co2_low = co2/s = co2/(co2_high/co2)
ind_toolow <- which(dat$co2_low == 0)
dat$co2_low[ind_toolow] <- (dat$co2[ind_toolow]^2)/dat$co2_high[ind_toolow]

# for that matter, get rid of any data points with central estimate of
# pCO2 <= 100
ind_remove <- which(dat$co2 <= 100)
dat <- dat[-ind_remove,]

# want to determine the two gamma parameters (shape, scale; mean=shape*scale)
# that give the median and 68% probability mass between upper/lower erorr bounds
# (dat$co2_low and dat$co2_high, with dat$co2 as the mean)

# define a function for fitting the parameters within a given tolerance
# underconstrained system, with 3 quantiles/mean to fit and 2 parameters
# Trying to fit all three gives exceedingly poor results, so only fit the two
# quantiles
# parameters[1]=shape
# parameters[2]=scale
sse_quantiles <- function(parameters, q16, q84) {
    q16.fit <- qlnorm(0.16, meanlog=parameters[1], sdlog=parameters[2])
    q84.fit <- qlnorm(0.84, meanlog=parameters[1], sdlog=parameters[2])
    sse <- sqrt((q16.fit-q16)^2 + (q84.fit-q84)^2)
    return(sse)
}

# preliminary testing to see how many iterations are needed to get convergence
# both shape and scale parameters (k and theta, with mean = k*theta, var = k*theta^2)
# are positive
bound.lower <- c(0, 0)
bound.upper <- c(1e5, 10)

# these settings yield maximum SSE of less than 0.01 ppmv (pCO2 quantiles)
niter.deoptim=1000        # number of iterations for DE optimization
NP.deoptim=50             # population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8             # as suggested by Storn et al (2006)
CR.deoptim=0.9            # as suggested by Storn et al (2006)

parameters.co2 <- vector('list',2)
names(parameters.co2) <- c('meanlog','sdlog')
for (p in 1:2) {parameters.co2[[p]] <- rep(NA,length(dat$co2))}
error.co2 <- rep(NA,length(dat$co2))
t1 <- proc.time()
pb <- txtProgressBar(min=0,max=length(dat$co2),initial=0,style=3)
for (i in 1:length(dat$co2)){
  outDEoptim <- DEoptim(sse_quantiles, bound.lower, bound.upper,
                      DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,
                      CR=CR.deoptim,trace=FALSE),
                      q16=dat$co2_low[i], q84=dat$co2_high[i])
  parameters.co2$meanlog[i] <- outDEoptim$optim$bestmem[1]
  parameters.co2$sdlog[i] <- outDEoptim$optim$bestmem[2]
  error.co2[i] <- sse_quantiles(outDEoptim$optim$bestmem, dat$co2_low[i], dat$co2_high[i])
  setTxtProgressBar(pb, i)
}
close(pb)
t2 <- proc.time()
plot(error.co2)

# check the worst error
i <- which.max(error.co2)
print(c(qlnorm(p=0.16, meanlog=parameters.co2$meanlog[i], sdlog=parameters.co2$sdlog[i]),qlnorm(p=0.84, meanlog=parameters.co2$meanlog[i], sdlog=parameters.co2$sdlog[i])))
print(dat[i,c('co2_low','co2','co2_high')])

# save workspace image
save.image(file = "fit_co2_data_logNormal.RData")

# write new CSV file, copy of Foster calib data set, but with the two log-normal
# parameters for each CO2 data point.
dat.co2 <- cbind(dat,parameters.co2$meanlog, parameters.co2$sdlog)
colnames(dat.co2) <- c(colnames(dat),'meanlog_co2','sdlog_co2')

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.out <- paste('../input_data/CO2_Proxy_Foster2017_calib_LN-co2_',today,'.csv',sep='')
write.csv(dat.co2, file=filename.out)

##==============================================================================
## End
##==============================================================================
