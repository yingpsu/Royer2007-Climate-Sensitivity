##==============================================================================
## plotting.R
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

library(ncdf4)

filename.mcmc <- 'geocarb_calibratedParameters_noPaleosols_03Jun2017.nc'
ncdata <- nc_open(filename.mcmc)
  parameters.nopaleosols = t(ncvar_get(ncdata, 'geocarb_parameters'))
  parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)

filename.mcmc <- 'geocarb_calibratedParameters_withPaleosols_03Jun2017.nc'
ncdata <- nc_open(filename.mcmc)
  parameters.paleosols = t(ncvar_get(ncdata, 'geocarb_parameters'))
  parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)



##==============================================================================


##==============================================================================
## End
##==============================================================================
