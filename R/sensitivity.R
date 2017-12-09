##==============================================================================
## sensitivity_geocarb.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())

co2_uncertainty_cutoff <- 20


if(Sys.info()['nodename']=='Tonys-MacBook-Pro.local') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
  .Ncore <- 0
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GECARB/R')
  .Ncore <- 15  # use multiple cores to process the many tide gauge stations
}



##==============================================================================
## Data
##=====

source('GEOCARB-2014_getData.R')

# remove the lowest [some number] co2 content data points (all paleosols, it turns out)
# (lowest ~40 are all from paleosols, actually)
#ind_co2_sort_all <- order(data_calib_all$co2)
#n_cutoff <- length(which(data_calib_all$co2 < quantile(data_calib_all$co2, 0.01)))
#data_calib_all <- data_calib_all[-ind_co2_sort_all[1:n_cutoff], ]

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , FALSE),
                        c("alkenones" , FALSE),
                        c("stomata"   , TRUE),
                        c("boron"     , FALSE),
                        c("liverworts", FALSE) )

ind_data    <- which(data_to_assim[2,]==TRUE)
n_data_sets <- length(ind_data)
ind_assim   <- vector("list",n_data_sets)
for (i in 1:n_data_sets) {
  ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
}

data_calib <- data_calib_all[unlist(ind_assim),]

# possible filtering out of some data points with too-narrow uncertainties in
# co2 (causing overconfidence in model simulations that match those data points
# well)
if(co2_uncertainty_cutoff > 0) {
  co2_halfwidth <- 0.5*(data_calib$co2_high - data_calib$co2_low)
  ind_filter <- which(co2_halfwidth < co2_uncertainty_cutoff)
  data_calib <- data_calib[-ind_filter,]
}

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
# the sensitivity analysis call below.
##==============================================================================


##==============================================================================
## Set up for optimization
##========================

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# need the likelihood function and prior distributions
source('GEOCARB-2014_calib_likelihood.R')
##==============================================================================


##==============================================================================
## Run the sensitivity analysis
##=============================

#install.packages('sensitivity')
library(sensitivity)
library(sn)

## Method of Morris
## Send in bounds of [0,1] for each parameter, and then scale within the
## `log_like_sensitivity` function to the parameters' distributions (using CDF)

mom <- morris(model=log_like_sensitivity, factors=parnames_calib, r=1000,
              binf=rep(0.01,length(parnames_calib)), bsup=rep(0.99,length(parnames_calib)),
              scale=FALSE, design=list(type="oat", levels=500, grid.jump=3),
              par_fixed=par_fixed0, parnames_calib=parnames_calib,
              parnames_fixed=parnames_fixed, age=age, ageN=ageN,
              ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
              ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
              input=input,
              data_calib=data_calib, ind_mod2obs=ind_mod2obs,
              ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
              iteration_threshold=iteration_threshold)

##==============================================================================








if(FALSE){



##==============================================================================
## Method of Morris

mom <- morris(model=neg_log_likelihood, factors=parnames, r = 20, binf=bound_lower, bsup=bound_upper,
              scale=TRUE, design=list(type="oat", levels=5, grid.jump=3),
              parnames=parnames, data_sealevel=data_sealevel,
              indices=indices, temperature_forcing=forcing)
##==============================================================================


##==============================================================================
## Sobol' Method

#install.packages('lhs')
library(lhs)

n_sample <- 1000
n_bootstrap <- 100

## Need preliminary function to map ranges from [0,1] and back
map_range <- function(X, lbin, ubin, lbout, ubout){
    Y <- lbout + (ubout-lbout)*( (X-lbin)/(ubin-lbin) )
    return(Y)
}

## Sample parameters (need 2 data frames)
parameters_lhs1 <- randomLHS(n_sample, n_parameters)
parameters_lhs2 <- randomLHS(n_sample, n_parameters)

## Sobol' wrapper - assumes uniform distributions on parameters
gmsl_sobol <- function(parameters_lhs, parnames, data_sealevel, indices,
                       temperature_forcing, bound_lower, bound_upper) {
  parameters <- sapply(1:length(parnames), function(pp) {
                  map_range(parameters_lhs[,pp], 0, 1, bound_lower[pp], bound_upper[pp])})
  nll_output <- neg_log_likelihood(parameters, parnames, data_sealevel, indices, temperature_forcing)
  nll_output_centered <- nll_output - mean(nll_output)
  return(nll_output_centered)
}

## Actually run the Sobol'
t.out <- system.time(s.out <- sobolSalt(model=gmsl_sobol,
                             parameters_lhs1,
                             parameters_lhs2,
                             scheme='B',
                             nboot=n_bootstrap))


##==============================================================================


##==============================================================================
## Latin Hypercube

todo
##==============================================================================

}

##==============================================================================
## End
##==============================================================================
