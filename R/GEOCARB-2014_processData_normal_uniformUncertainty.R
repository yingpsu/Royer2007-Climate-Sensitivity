##==============================================================================
## GEOCARB-2014_processData_normal_uniformUncertainty.R
##
## Read CO2 proxy data, with age/amount uncertainties.
##
## Fit normal distributions at each time point, assuming the +sigma1/-sigma2
## range contains 0.68 probability mass. Want to fit the 16%, 84% quantiles and
## the mean.
##
## Write out a file with the normal distribution parameters on it, to use
## in the calibration.
##
## --> this version sets up experiment using a single uncertainty for all
##     of the data points, to assess the extent to which the GEOCARB low-CO2
##     bias is attributable to higher uncertainties for higher-CO2 data points
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

# col 1 is age (million years ago), col 2 is CO2 (units?), col 3 is lower error
# bound, col 4 is upper error bound
dat <- read.csv('../input_data/CO2_Proxy_Foster2017_calib.csv', fill=TRUE, header=TRUE)

parameters.co2 <- vector('list',2)
names(parameters.co2) <- c('mu','sigma')

# take normal parameter mu to be the central estimate
parameters.co2$mu <- dat$co2
parameters.co2$sigma <- 0.5*( (dat$co2_high-dat$co2) + (dat$co2-dat$co2_low))

parameters.co2$sigma <- rep(mean(parameters.co2$sigma), length(parameters.co2$sigma))

# write new CSV file, copy of Foster calib data set, but with the two normal
# parameters for each CO2 data point.
dat.co2 <- cbind(dat,parameters.co2$mu, parameters.co2$sigma)
colnames(dat.co2) <- c(colnames(dat),'mu_co2','sigma_co2')

# clean out the rows that have co2 <= 0, or co2_low=co2_high=0. there are not useful.
ind_co2 <- which(dat.co2$co2 <= 0)
ind_co2_low <- which(dat.co2$co2_low==0)
ind_co2_high <- which(dat.co2$co2_high==0)
irem <- unique(c( intersect(ind_co2_low,ind_co2_high),ind_co2))
dat.co2 <- dat.co2[-irem,]

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.out <- paste('../input_data/CO2_Proxy_Foster2017_calib_NM-co2-unifUnc_',today,'.csv',sep='')
write.csv(dat.co2, file=filename.out)

##==============================================================================
## End
##==============================================================================
