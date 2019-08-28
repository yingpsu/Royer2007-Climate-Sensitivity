##==============================================================================
## calib_driver_nm-mixture.R
##
## Read CO2 proxy data, with fitted parameters for skew-normal kernels for each
## data point. Set which data sets you intend to calibrate using.
## Version for running as script on HPC.
##
## Uses the fit_likelihood_surface_mixture.R script and the `dist` variable to
## fit a mixture model based on kernels of the form `dist`.
##
## Values are currently set for a calibration that will be feasible on a modern
## laptop computer (< 1 hr).
## File names will need to be modified to fit what you have named your files.
##
## This version uses the normal kernels to fit a Gaussian mixture model for
## each time slice in the likelihood function (as opposed to the skew-normal
## mixture model of the main results).
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

print(paste("START AT",Sys.time()))

rm(list=ls())

setwd('~/work/codes/GEOCARB/R')

niter_mcmc000 <- 1e4   # number of MCMC iterations per node (Markov chain length)
n_node000 <- 1        # number of CPUs to use
appen <- 'mix' # 'mix' for main results; 'PR2011' for supplemental experiment; 'sens' for only calibrating sensitive parameters (after plotting_sobol.R)
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")

# Distribution fit to each data point in the processing step
dist <- 'nm-mix'  # normal (use this to reproduce supplementary experiment results)
#dist <- 'sn-mix'  # skew-normal mixture model
appen2 <- dist

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- FALSE
DO_PARAM_INIT <- FALSE # do initialization of parameters & covariance matrix from previous calibration?
USE_LENTON_FSR <- FALSE
USE_ROYER_FSR <- TRUE

filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.covarinit <- "../output/covar_init_sn-mix_12Aug2019.rds"
filename.paraminit <- "../output/param_init_sn-mix_12Aug2019.rds"

library(adaptMCMC)
library(sn)
library(invgamma)

##==============================================================================
## Model parameters and setup
##===========================

source('parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')

#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `stdev` uncertainty statistical parameter
##==============================================================================


##==============================================================================
## Initialization of  parameters and transition matrix
##====================================================
# quick fix to initialize standard deviation
if(!is.na(match('stdev',parnames_calib))) {
  par_calib0[match('stdev',parnames_calib)] <- 450
}

if(DO_PARAM_INIT) {
  step_mcmc <- readRDS(filename.covarinit)
  par_calib0 <- readRDS(filename.paraminit) # this will overwrite setting stdev above
}
##==============================================================================


##==============================================================================
## Data
##=====

source('fit_likelihood_surface_mixture.R')
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
stopadapt_mcmc <- round(niter_mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)

##==============================================================================
## Actually run the calibration

if(n_node000==1) {
  tbeg <- proc.time()
  amcmc_out1 <- MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                  scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=max(5000,round(0.05*niter_mcmc)),
                  par_fixed=par_fixed0, parnames_calib=parnames_calib,
                  parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                  ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                  ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                  input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                  data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                  ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                  iteration_threshold=iteration_threshold,
                  loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                  do_sample_tvq=TRUE, par_time_center=par_time_center, par_time_stdev=par_time_stdev,
                  upper_bound_co2=.upper_bound_co2, lower_bound_co2=.lower_bound_co2)
  tend <- proc.time()
  chain1 = amcmc_out1$samples
} else if(n_node000 > 1) {
  tbeg <- proc.time()
  amcmc_par1 <- MCMC.parallel(log_post, n=niter_mcmc, init=par_calib0, n.chain=n_node000, n.cpu=n_node000,
        					dyn.libs=c('../fortran/run_geocarb_unc.so'),
                  packages=c('sn','invgamma'),
        					adapt=TRUE, list=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
        					gamma=gamma_mcmc, n.start=max(5000,round(0.05*niter_mcmc)),
                  par_fixed=par_fixed0, parnames_calib=parnames_calib,
                  parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                  ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                  ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                  input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                  data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                  ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                  iteration_threshold=iteration_threshold,
                  loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                  do_sample_tvq=TRUE, par_time_center=par_time_center, par_time_stdev=par_time_stdev,
                  upper_bound_co2=.upper_bound_co2, lower_bound_co2=.lower_bound_co2)
  tend <- proc.time()
}
print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

## write an RData file with the parameter results
filename.mcmc = paste(output_dir,'geocarb_mcmcoutput_',appen,'_',today,appen2,'.RData',sep="")
if(n_node000 == 1) {
  save(amcmc_out1, file=filename.mcmc)
} else {
  save(amcmc_par1, file=filename.mcmc)
}

## Extend and run more MCMC samples?
if(FALSE){
niter_extend <- 5e6
tbeg <- proc.time()
amcmc_extend1 <- MCMC.add.samples(amcmc_out1, niter_extend,
                                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                                parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                                input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                                data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                                iteration_threshold=iteration_threshold,
                                loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                                do_sample_tvq=TRUE, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
tend <- proc.time()
chain1 <- amcmc_extend1$samples
amcmc_out1 <- amcmc_extend1
}
##==============================================================================


##==============================================================================
## Convergence diagnostics and autocorrelation accounted for in the
## burnin_thinning.R script. So we're done!

print(paste("END AT",Sys.time()))

##==============================================================================


##==============================================================================
## End
##==============================================================================
