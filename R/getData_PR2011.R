##==============================================================================
## getData_PR2011.R
##
## Read CO2 proxy data, including skew-normal/normal fits for CO2 uncertainty
## distributions (previously fit using processData_[something].R).
##
## Assumes that `filename.data` will have been set in the calling routine, and
## customized to read the Park and Royer 2011 data set.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

data_calib_all <- read.csv(filename.data, fill=TRUE, header=TRUE)

ind_assim <- 1:nrow(data_calib_all)
data_calib <- data_calib_all[unlist(ind_assim),]

# assumption of steady state in-between model time steps permits figuring out
# which model time steps each data point should be compared against in advance.
# doing this each calibration iteration would be outrageous!
# This assumes the model time step is 10 million years, seq(570,0,by=-10). The
# model will choke later (in calibration) if this is not consistent with what is
# set within the actual GEOCARB physical model.
age_tmp <- seq(570,0,by=-10)
ttmp <- 10*ceiling(data_calib$age/10)
ind_mod2obs <- rep(NA,nrow(data_calib))
for (i in 1:length(ind_mod2obs)){
  ind_mod2obs[i] <- which(age_tmp==ttmp[i])
}

##==============================================================================
## End
##==============================================================================
