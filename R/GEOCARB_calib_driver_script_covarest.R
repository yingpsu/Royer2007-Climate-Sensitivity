##==============================================================================
## GEOCARB-2014_calib_driver_script_covarest.R
##
## Read CO2 proxy data. Set which data sets you intend to calibrate using.
## Version for running as script on HPC.
##
## Will run subsets of the calibration parameters to obtain a better starting
## estimate of the covariance matrix for the whole bunch, and better initial
## parameter estimates.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

niter_mcmc000 <- 4e4   # number of MCMC iterations to use for each subsample
n_node000 <- 1         # number of CPUs to use
setwd('/home/scrim/axw322/codes/GEOCARB/R')
#setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
appen <- 'covarest'
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
l_write_rdata  <- FALSE
l_write_netcdf <- FALSE
co2_uncertainty_cutoff <- 20

library(sn)
library(adaptMCMC)
library(ncdf4)
library(magic)

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
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

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
# set to +65%, - 30% uncertain range around the central estimate
if(co2_uncertainty_cutoff > 0) {
  co2_halfwidth <- 0.5*(data_calib$co2_high - data_calib$co2_low)
  ind_filter <- which(co2_halfwidth < co2_uncertainty_cutoff)
  ind_remove <- NULL
  for (ii in ind_filter) {
    range_original <- data_calib[ii,'co2_high']-data_calib[ii,'co2_low']
    range_updated  <- data_calib[ii,'co2']*0.95
    if (range_updated > range_original) {
      # update to the wider uncertain range if +65/-30% is wider
      data_calib[ii,'co2_high'] <- data_calib[ii,'co2']*1.65
      data_calib[ii,'co2_low']  <- data_calib[ii,'co2']*0.70
    } else {
      # otherwise, remove
      ind_remove <- c(ind_remove, ii)
    }
  }
  data_calib <- data_calib[-ind_filter,]
  ##data_calib <- data_calib[-ind_remove,]
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

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# need the likelihood function and prior distributions
source('GEOCARB-2014_calib_likelihood.R')

# Read parameter information, set up the calibration parameters
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all.csv'
calib_all <- read.csv(filename.calibinput)
source('GEOCARB-2014_parameterSetup.R')
parnames_calib_all <- parnames_calib
n_params <- length(parnames_calib_all)
##==============================================================================


##==============================================================================
## Divide into parameter subsets, and calibrate
##=============================================

# break up into subsets.
# first subset is the expert set: {ACT, LIFE, GYM, FERT, GLAC, deltaT2X}
# other subsets are of size n_subset (50 parameters left, so do 10 subsets of 5)
p_expert <- c('ACT','LIFE','GYM','FERT','deltaT2X','GLAC')
n_subset <- 5
n_sets <- 1+ceiling((n_params-length(p_expert))/n_subset)  # these codes won't work if not divisible
covars <- vector('list', n_sets)
params <- rep(0, n_params)
amcmc <- vector('list', n_sets)

# set up which indices are in which subset
ind_subset <- vector('list', n_sets)
all_indices <- 1:n_params # and remove the ones we have assigned to a subset
ind_subset[[1]] <- match(p_expert, parnames_calib_all)
all_indices <- all_indices[-match(ind_subset[[1]], all_indices)]
for (k in 2:n_sets) {
  ind_subset[[k]] <- all_indices[1:n_subset]
  all_indices <- all_indices[-match(ind_subset[[k]], all_indices)]
}

for (k in 1:n_sets) {

  ##============================================
  ## Set up the parameter subset for calibration
  ##============================================

  # start by setting all of the calibration parameter flags to 0
  calib_sub <- calib_all
  calib_sub$calib <- 0

  # then set the ones we want to calibrate to 1, and write the temporary file
  for (j in 1:length(ind_subset[[k]])) {
    row <- which(calib_sub$parameter==parnames_calib_all[ind_subset[[k]][j]])
    calib_sub$calib[row] <- 1
  }
  write.csv(x=calib_sub, file='../input_data/GEOCARB_input_summaries_calib_subsample.csv', row.names=FALSE)

  # read those parameters and set up to calibrate em
  filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_subsample.csv'
  source('GEOCARB-2014_parameterSetup.R')

  # NB:  the parameters in each subset need to be in order (increasing index)

  # proceed with calibration...

  ##==========================================
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

  ##====================
  ## Run the calibration
  ##====================

  # interpolate between lots of parameters and one parameter.
  # this functional form yields an acceptance rate of about 25% for as few as 10
  # parameters, 44% for a single parameter (or Metropolis-within-Gibbs sampler),
  # and 0.234 for infinite number of parameters, using accept_mcmc_few=0.44 and
  # accept_mcmc_many=0.234.
  accept_mcmc_few <- 0.44         # optimal for only one parameter
  accept_mcmc_many <- 0.234       # optimal for many parameters
  accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_calib)
  niter_mcmc <- niter_mcmc000
  gamma_mcmc <- 0.5
  stopadapt_mcmc <- round(niter_mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)

  ## Actually run the calibration

  tbeg=proc.time()
  amcmc[[k]] = MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                  scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=max(1000,round(0.05*niter_mcmc)),
                  par_fixed=par_fixed0, parnames_calib=parnames_calib,
                  parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                  ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                  ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                  input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                  data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                  ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                  iteration_threshold=iteration_threshold)
  tend=proc.time()
  chain1 = amcmc[[k]]$samples
  print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

  ## Save the results
  params[ind_subset[[k]]] <- chain1[niter_mcmc,]
  covars[[k]] <- amcmc[[k]]$cov.jump
  rownames(covars[[k]]) <- colnames(covars[[k]]) <- parnames_calib

}
##==============================================================================


##==============================================================================
## Put back together into one covariance matrix estimate
##======================================================

# this assumes independence, which is of course not true.  but it gives a better
# initial estimate of the transition matrix than the `step` estimates

# need to put back in the proper order, or will indexing here take care of that?

covar_full <- adiag(covars[[1]], covars[[2]])
if (n_sets > 2) {
  for (k in 3:n_sets) {
    covar_full <- adiag(covar_full, covars[[k]])
  }
}
# rearrange back into original order
covar <- covar_full[parnames_calib_all, parnames_calib_all]

names(params) <- parnames_calib_all

# Save results to use later for better initial estimates
save(list=c('amcmc'), file='markov_chains.RData')
save(list=c('params','covar_full'), file='initial_estimates.RData')
##==============================================================================


##==============================================================================
## End
##==============================================================================
