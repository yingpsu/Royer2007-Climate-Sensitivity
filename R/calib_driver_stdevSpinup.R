##==============================================================================
## calib_driver_stdevSpinup.R
##
## Preliminary 1D calibration to spin up the uncertainty parameter and obtain
## good starting values and step sizes for it.
##
## Using the PR2011 data set and unimodal, with the old Royer et al 2014 fSR
## time series.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

print(paste("START AT",Sys.time()))

rm(list=ls())

setwd('~/work/codes/GEOCARB/R')

param_choice <- 'stdevOnly'   # Calibrate all 69 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'PR2011'    # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
lhood_choice <- 'unimodal'  # Mixture model ("mixture") or unimodal ("unimodal")?
dist <- 'sn'               # kernel choice for each data point (sn (skew-normal), ln (log-normal), nm (normal))

today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.mcmc = paste('../output/geocarb_mcmcoutput_stdevSpinup_',today,'.RData',sep="")

niter_mcmc000 <- 5000   # number of MCMC iterations per node (Markov chain length)
n_node000 <- 1        # number of CPUs to use

# If using the Foster et al 2017 data set, which proxy sets to assimilate?
#   (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

USE_LENTON_FSR <- FALSE
USE_ROYER_FSR <- FALSE

# do initialization of parameters & covariance matrix from previous calibration?
DO_PARAM_INIT <- TRUE   # if true, set the two file names below
filename.covarinit <- "../output/covar_init_sn-mix_12Aug2019.rds"
filename.paraminit <- "../output/param_init_sn-mix_12Aug2019.rds"

library(adaptMCMC)
library(sn)
library(invgamma)

##==============================================================================
## Model parameters and setup
##===========================

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

# sampling the time-varying arrays
DO_SAMPLE_TVQ <- FALSE

# Now, set up the parameters as desired. Possibly the same.
filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',param_choice,'.csv', sep='')
source('parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')

#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `stdev` uncertainty statistical parameter
##==============================================================================


##==============================================================================
## Initialization of  parameters and transition matrix
##====================================================

par_calib0 <- 120
step_mcmc <- 100
##==============================================================================


##==============================================================================
## Data
##=====

source('fit_likelihood_surface.R')
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

rm(list=c('bound_lower','bound_upper','bounds'))

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the calibration MCMC call below.
##==============================================================================


##==============================================================================
## Run the calibration
##====================

# need the likelihood function and prior distributions
source('calib_likelihood_unc.R')

# set up and run the actual calibration
# interpolate between lots of parameters and one parameter.
# this functional form yields an acceptance rate of about 25% for as few as 10
# parameters, 44% for a single parameter (or Metropolis-within-Gibbs sampler),
# and 0.234 for infinite number of parameters, using accept_mcmc_few=0.44 and
# accept_mcmc_many=0.234.
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames_calib)
niter_mcmc <- niter_mcmc000
gamma_mcmc <- 0.66
startadapt_mcmc <- 1000

##==============================================================================
## Actually run the calibration

print(paste("data =",data_choice," | params =",param_choice," | lhood =",lhood_choice))
tbeg <- proc.time()
amcmc_out1 <- MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                  scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                  par_fixed=par_fixed0, parnames_calib=parnames_calib,
                  parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                  ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                  ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                  input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                  data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                  ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                  iteration_threshold=iteration_threshold,
                  loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                  do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev,
                  upper_bound_co2=.upper_bound_co2, lower_bound_co2=.lower_bound_co2)
tend <- proc.time()
chain1 = amcmc_out1$samples

print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

## write an RData file with the parameter results
save(amcmc_out1, file=filename.mcmc)
##==============================================================================


##==============================================================================
## Convergence diagnostics and autocorrelation accounted for in the
## burnin_thinning.R script. So we're done!

print(paste("END AT",Sys.time()))

##==============================================================================


##==============================================================================
## End
##==============================================================================
