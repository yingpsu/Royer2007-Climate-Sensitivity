##==============================================================================
## supplemental_experiment_PR2011_unimodalfit.R
##
## Run a supplemental calibration experiment where only the 6 parameters
## from Park and Royer (2011) are varied: ACT, GYM, FERT, LIFE, deltaT2X, GLAC.
##
## Also, fitting only a single Gaussian (truncated in likelihood function)
## distribution for each time slice in processing step. To assess the impacts
## of our improved error model.
##
## The processing and model ensemble simulations are all done in this script
## too because it's a small experiment and can be run on a modern laptop.
##
## Yields output:
##   geocarb_mcmcoutput_PR2011_[DATESTAMP]unimodal.RData (only the MCMC output)
##   GEOCARB_MCMC_PR2011_[DATESTAMP]unimodal.RData (full workspace from MCMC)
##   geocarb_calibratedParameters_PR2011_[DATESTAMP]unimodal.RData (just the parameters after post-processing)
##   analysis_PR2011_[DATESTAMP]unimodal.RData (output + analysis)
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

print(paste("START AT",Sys.time()))

rm(list=ls())

setwd('~/codes/GEOCARB/R') # set up for Tony's machine; yours is probably different... :)

niter_mcmc000 <- 2e5   # number of MCMC iterations per node (Markov chain length)
n_node000 <- 1        # number of CPUs to use
appen <- 'PR2011' # 'unc' for main results; 'PR2011' for supplemental experiment
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")

# Distribution fit to each data point in the processing step
#dist <- 'ga'  # gamma
#dist <- 'be'  # beta
#dist <- 'ln'  # log-normal
#dist <- 'sn'  # skew-normal (use this to reproduce main results)
#dist <- 'nm'  # normal (use this to reproduce supplementary experiment results)
#dist <- 'sn-100min'  # skew-normal (use this to reproduce supplementary experiment results)
#dist <- 'sn-mmrem'  # skew-normal (use this to reproduce supplementary experiment results)
#dist <- 'nm-unifUnc' # normal (but with all data points assigned the same uncertainty)
dist <- 'uni' # unimodal experiment
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

DO_SAMPLE_TVQ <- FALSE  # sample time series uncertainty by CDF parameters?
DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- FALSE
DO_PARAM_INIT <- FALSE # do initialization of parameters & covariance matrix from previous calibration?
USE_LENTON_FSR <- FALSE
USE_DT2019_FSR <- TRUE

filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.covarinit <- "../output/covar_init_unc-sd10_22Jun2019.rds"
filename.paraminit <- "../output/param_init_unc-sd10_22Jun2019.rds"

library(adaptMCMC)
library(sn)
library(invgamma)

##==============================================================================
## Model parameters and setup
##===========================

if(DO_SAMPLE_TVQ) {
  source('GEOCARB-2014_parameterSetup_tvq.R')
  source('model_forMCMC_tvq.R')
} else {
  source('GEOCARB-2014_parameterSetup.R')
  source('model_forMCMC.R')
}

#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `var` uncertainty statistical parameter
##==============================================================================


##==============================================================================
## Initialization of  parameters and transition matrix
##====================================================
# quick fix to initialize standard deviation
if(!is.na(match('stdev',parnames_calib))) {
  par_calib0[match('stdev',parnames_calib)] <- 450
  #par_calib0[match('stdev',parnames_calib)] <- rinvgamma(shape=input[input$parameter=='var', 'mean'], rate=input[input$parameter=='var', 'two_sigma'], n=1)
}

if(DO_PARAM_INIT) {
  step_mcmc <- readRDS(filename.covarinit)
  par_calib0 <- readRDS(filename.paraminit) # this will overwrite setting stdev above
}
##==============================================================================


##==============================================================================
## Data
##=====

source('GEOCARB_fit_likelihood_surface_unimodal.R')
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
source('GEOCARB-2014_calib_likelihood_unc.R')

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
                  upper_bound_co2=.upper_bound_co2, lower_bound_co2=.lower_bound_co2)
  tend <- proc.time()
}

print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

# save
filename.out <- paste(output_dir,'GEOCARB_MCMC_',appen,'_',today,appen2,'.RData', sep='')
if(DO_WRITE_RDATA) {save.image(file=filename.out)}

## write an RData file with the parameter results
filename.mcmc = paste(output_dir,'geocarb_mcmcoutput_',appen,'_',today,appen2,'.RData',sep="")
if(n_node000 == 1) {
  save(amcmc_out1, file=filename.mcmc)
} else {
  save(amcmc_par1, file=filename.mcmc)
}

##==============================================================================


##==============================================================================
## Convergence diagnostics
##========================

## Gelman and Rubin diagnostics - determine and chop off for burn-in
initial <- 0.5*niter_mcmc
niter.test <- initial
gr.test <- rep(0, length(niter.test))

if(n_node000 == 1) {
    # don't do GR stats, just cut off first half of chains
    print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
} else if(n_node000 > 1) {
    # this case is FAR more fun
    # accumulate the names of the soon-to-be mcmc objects
    string.mcmc.list <- 'mcmc1'
    for (m in 2:n_node000) {
        string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
    }
    for (i in 1:length(niter.test)) {
        for (m in 1:n_node000) {
            # convert each of the chains into mcmc object
            eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_par1[[m]]$samples[(initial+1):niter_mcmc,])', sep='')))
        }
        eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

        gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
    }
} else {print('error - n_node000 < 1 makes no sense')}

##==============================================================================




##==============================================================================
## Chop off burn-in
##==============================================================================

ifirst <- NA
if(n_node000==1) {
  ifirst <- round(0.5*niter_mcmc)
} else {
  gr.max <- 1.1
  for (i in seq(from=length(niter.test), to=1, by=-1)) {
    if( all(gr.test[i:length(gr.test)] < gr.max) ) {ifirst <- niter.test[i]}
  }
}


chains_burned <- NA
if(n_node000 > 1) {
  chains_burned <- vector('list', n_node000)
  for (m in 1:n_node000) {
    chains_burned[[m]] <- amcmc_par1[[m]]$samples[(ifirst+1):niter_mcmc,]
  }
} else {
  chains_burned <- amcmc_out1$samples[(ifirst+1):niter_mcmc,]
}

##==============================================================================
## possible thinning?
##==============================================================================



lmax <- 200
cmax <- 0.05
maxlag <- 0

for (m in 1:n_node000) {
    for (p in 1:length(parnames_calib)) {
        acf_tmp <- acf(chains_burned[[m]][,p], lag.max=lmax, plot=FALSE)
        idx_low <- which(acf_tmp$acf < cmax)
        while (length(idx_low)==0) {
          lmax <- lmax + 100
          acf_tmp <- acf(chains_burned[[m]][,p], lag.max=lmax, plot=FALSE)
          idx_low <- which(acf_tmp$acf < cmax)
        }
        new <- acf_tmp$lag[idx_low[1]]
        if (maxlag < new) {
            print(paste(m,p,"Updating maxlag to",new))
            maxlag <- new
        }
    }
}

chains_burned_thinned <- chains_burned  # initialize
if(n_node000 > 1) {
  for (m in 1:n_node000) {
    chains_burned_thinned[[m]] <- chains_burned[[m]][seq(from=1, to=nrow(chains_burned[[m]]), by=maxlag),]
  }
} else {
  chains_burned_thinned <- chains_burned[seq(from=1, to=nrow(chains_burned), by=maxlag),]
}


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# samples from the posterior. Only saving the transition covariance matrix for
# one of the chains (if in parallel).

if(n_node000==1) {
  parameters.posterior <- chains_burned_thinned
  covjump.posterior <- amcmc_out1$cov.jump
} else {
  parameters.posterior <- chains_burned_thinned[[1]]
  covjump.posterior <- amcmc_par1[[1]]$cov.jump
  for (m in 2:n_node000) {
    parameters.posterior <- rbind(parameters.posterior, chains_burned_thinned[[m]])
  }
}


##==============================================================================
## Write calibrated parameters output file
##========================================

filename.parameters = paste(output_dir,'geocarb_calibratedParameters_',appen,'_',today,appen2,'.RData',sep="")
save(parameters.posterior, parnames_calib, covjump.posterior, file=filename.parameters)

##==============================================================================



##==============================================================================
## Run the model ensemble
##=======================

# run the ensemble
model_pr2011 <- sapply(X=1:nrow(parameters.posterior),
              FUN=function(k){model_forMCMC(par_calib=parameters.posterior[k,],
                                            par_fixed=par_fixed0,
                                            parnames_calib=parnames_calib,
                                            parnames_fixed=parnames_fixed,
                                            parnames_time=parnames_time,
                                            age=age,
                                            ageN=ageN,
                                            ind_const_calib=ind_const_calib,
                                            ind_time_calib=ind_time_calib,
                                            ind_const_fixed=ind_const_fixed,
                                            ind_time_fixed=ind_time_fixed,
                                            ind_expected_time=ind_expected_time,
                                            ind_expected_const=ind_expected_const,
                                            iteration_threshold=iteration_threshold,
                                            do_sample_tvq=DO_SAMPLE_TVQ,
                                            par_time_center=par_time_center,
                                            par_time_stdev=par_time_stdev)[,'co2']})
n_time <- nrow(model_pr2011)

# get 5-95% range and median  are cols 1-3; max-post will be 4
quantiles_i_want <- c(0,0.005,.025,.05,.5,.95,.975,0.995,1)
model_quantiles <- mat.or.vec(nr=n_time, nc=(length(quantiles_i_want)+1))
colnames(model_quantiles) <- c('q000','q005','q025','q05','q50','q95','q975','q995','q100','maxpost')
for (t in 1:n_time) {
    model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_pr2011[t,], quantiles_i_want)
}

# get posterior scores

lpost_out <- sapply(X=1:nrow(parameters.posterior),
              FUN=function(k){log_post(par_calib=parameters.posterior[k,],
                                      par_fixed=par_fixed0,
                                      parnames_calib=parnames_calib,
                                      parnames_fixed=parnames_fixed,
                                      parnames_time=parnames_time,
                                      age=age,
                                      ageN=ageN,
                                      ind_const_calib=ind_const_calib,
                                      ind_time_calib=ind_time_calib,
                                      ind_const_fixed=ind_const_fixed,
                                      ind_time_fixed=ind_time_fixed,
                                      input=input,
                                      time_arrays=time_arrays,
                                      bounds_calib=bounds_calib,
                                      data_calib=data_calib,
                                      ind_mod2obs=ind_mod2obs,
                                      ind_expected_time=ind_expected_time,
                                      ind_expected_const=ind_expected_const,
                                      iteration_threshold=iteration_threshold,
                                      n_shard=1,
                                      loglikelihood_smoothed=loglikelihood_smoothed,
                                      likelihood_fit=likelihood_fit,
                                      idx_data=idx_data,
                                      do_sample_tvq=DO_SAMPLE_TVQ,
                                      par_time_center=par_time_center,
                                      par_time_stdev=par_time_stdev)})

model_quantiles[,'maxpost'] <- model_pr2011[,which.max(lpost_out)]

parnames_pr2011uni <- parnames_calib
parameters_posterior_pr2011uni <- parameters.posterior
model_quantiles_pr2011uni <- model_quantiles

## save everything you might need in a special file
filename.analysis = paste(output_dir,'analysis_',appen,'_',today,'unimodal.RData',sep="")
save(parameters_posterior_pr2011uni, parnames_pr2011uni, model_quantiles_pr2011uni, file=filename.analysis)

## and save *everything* in the analysis file
if(DO_WRITE_RDATA) {save.image(file=filename.out)}

##==============================================================================


##==============================================================================
## End
##==============================================================================
