##==============================================================================
## GEOCARB-2014_calib_driver.R
##
## Read CO2 proxy data. Set which data sets you intend to calibrate using.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================


##==============================================================================
## Data
##=====

# Read proxy data
source('getData.R')

# Which proxy sets to use?
#TODO (once Ying puts another column to designate different outputs)
##==============================================================================


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
source('GEOCARB-2014_parameterSetup.R')
##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distributions
names <- as.character(input$parameter)
bound_lower <- rep(NA, length(names))
bound_upper <- rep(NA, length(names))

ind_neg_inf <- which(input[,'lower_limit']=='_inf')
bound_lower[ind_neg_inf] <- -Inf
bound_lower[setdiff(1:length(names), ind_neg_inf)] <- as.numeric(as.character(input$lower_limit[setdiff(1:length(names), ind_neg_inf)]))
bound_upper <- input$upper_limit

bounds <- cbind(bound_lower, bound_upper)
rownames(bounds) <- as.character(input$parameter)

# only actually need the calibration parameters' bounds, so reformat the bounds
# array to match the vector of calibration parameters
bounds_calib <- mat.or.vec(nr=length(parnames_calib), nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:length(parnames_calib)) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the calibration MCMC call below.
##==============================================================================


##==============================================================================
## Run the calibration
##====================


##==============================================================================


##==============================================================================
## Convergence diagnostics
##========================

#TODO
# Gelman and Rubin
# Heidelberger and Welch
# visual inspection
##==============================================================================


##==============================================================================
## Write calibrated parameters output file
##========================================


##==============================================================================


##==============================================================================
## End
##==============================================================================
