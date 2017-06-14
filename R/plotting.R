##==============================================================================
## plotting.R
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

library(ncdf4)

filename.mcmc <- 'geocarb_calibratedParameters_allData-allConst_12Jun2017.nc'
ncdata <- nc_open(filename.mcmc)
  parameters = t(ncvar_get(ncdata, 'geocarb_parameters'))
  parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)

filename.mcmc <- 'geocarb_calibratedParameters_withPaleosols_03Jun2017.nc'
ncdata <- nc_open(filename.mcmc)
  parameters.paleosols = t(ncvar_get(ncdata, 'geocarb_parameters'))
  parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)


par(mfrow=c(7,8))
for (p in 1:length(parnames)) {plot(parameters[,p], type='l', ylab=parnames[p])}

par(mfrow=c(1,1))
hist(parameters[round(nrow(parameters)*0.5):nrow(parameters), match('deltaT2X',parnames)],
     freq=FALSE, xlab='deltaT2X [deg C]', main='')
##==============================================================================


##==============================================================================
## End
##==============================================================================
