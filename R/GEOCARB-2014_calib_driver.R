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

# Read parameter information
input.summary <- read.csv("GEOCARB_input_summaries_calib.csv")

# Which model parameters to calibrate?
#TODO
# read column "calib" in GEOCARB_input_summaries_calib.csv; a 0 implies to keep
# the parameter/forcing fixed at the central estimate on the file, and a 1
# implies to calibrate the parameter

##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distributions
#TODO

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
