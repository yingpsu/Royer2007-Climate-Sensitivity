##==============================================================================
## GEOCARB-2014_calib_driver_split.R
##
## Read CO2 proxy data. Set which data sets you intend to calibrate using.
## Version for running as script on HPC.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')

niter_mcmc000 <- 2e6   # number of MCMC iterations per node (Markov chain length)
n_node000 <- 1        # number of CPUs to use
appen <- 'tvq_split'
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")

# Distribution fit to each data point in the processing step
#dist <- 'ga'  # gamma
#dist <- 'be'  # beta
#dist <- 'ln'  # log-normal
dist <- 'sn'  # skew-normal (use this to reproduce main results)
#dist <- 'nm'  # normal (use this to reproduce supplementary experiment results)
#dist <- 'sn-100min'  # skew-normal (use this to reproduce supplementary experiment results)
#dist <- 'sn-mmrem'  # skew-normal (use this to reproduce supplementary experiment results)
appen2 <- paste(dist,'-split-deg',sep='')

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

DO_SAMPLE_TVQ <- TRUE  # sample time series uncertainty by CDF parameters?
DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- TRUE
USE_LENTON_FSR <- TRUE

filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_unc.csv'
#filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')

library(adaptMCMC)
library(ncdf4)
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
source('run_geocarbF_unc.R')
##==============================================================================


# quick fix to initialize variance
par_calib0[match('var',parnames_calib)] <- rinvgamma(shape=input[input$parameter=='var', 'mean'], rate=input[input$parameter=='var', 'two_sigma'], n=1)


##==============================================================================
## Data
##=====

source('GEOCARB_fit_likelihood_surface_split.R')
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
source('GEOCARB-2014_calib_likelihood.R')

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

if(DO_SAMPLE_TVQ) {
  if(n_node000==1) {
    amcmc_out <- vector('list', n_shard)
    tbeg <- proc.time()
    for (ss in 1:n_shard) {
      amcmc_out[[ss]] <- MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                    scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=max(5000,round(0.05*niter_mcmc)),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                    ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                    ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                    input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                    data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                    ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                    iteration_threshold=iteration_threshold,
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit[[ss]], idx_data=idx_data[[ss]],
                    do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
    }
    tend <- proc.time()
    #chain1 = amcmc_out1$samples
  } else if(n_node000 > 1) {
    amcmc_par <- vector('list', n_shard)
    tbeg <- proc.time()
    for (ss in 1:n_shard) {
      amcmc_par[[ss]] <- MCMC.parallel(log_post, n=niter_mcmc, init=par_calib0, n.chain=n_node000, n.cpu=n_node000,
          					dyn.libs=c('../fortran/run_geocarb.so'),
                    packages=c('sn'),
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
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit[[ss]], idx_data=idx_data[[ss]],
                    do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
    }
    tend <- proc.time()
  }
} else {
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
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data)
    tend <- proc.time()
    chain1 = amcmc_out1$samples
  } else if(n_node000 > 1) {
    tbeg <- proc.time()
    amcmc_par1 <- MCMC.parallel(log_post, n=niter_mcmc, init=par_calib0, n.chain=n_node000, n.cpu=n_node000,
          					dyn.libs=c('../fortran/run_geocarb.so'),
                    packages=c('sn'),
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
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data)
    tend <- proc.time()
  }
}
print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

# save
if(DO_WRITE_RDATA) {save.image(file=paste(output_dir,'GEOCARB_MCMC_',appen,'_',today,appen2,'.RData', sep=''))}

## Plot a little?
#for (ss in 1:3) {for (mm in 1:4) {plot(amcmc_par[[ss]][[mm]]$samples[,10], type='l')}}

## Extend an MCMC chain?
## Extend and run more MCMC samples?
if(FALSE){
niter_extend <- 1e4
tbeg=proc.time()
amcmc_extend1 = MCMC.add.samples(amcmc_out1, niter_extend,
                                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                                parnames_fixed=parnames_fixed, parnames_time=parnames_time, age=age, ageN=ageN,
                                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                                input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                                data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                                iteration_threshold=iteration_threshold,
                                loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                                do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
tend=proc.time()
chain1 = amcmc_extend1$samples
}
##==============================================================================


##==============================================================================
## Convergence diagnostics
##========================

## Gelman and Rubin diagnostics - determine and chop off for burn-in
initial <- 0.5*niter_mcmc
increment <- round(0.05*niter_mcmc)
niter.test <- seq(from=(round(0.5*niter_mcmc)+increment), to=niter_mcmc, by=increment)
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
      eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_par1[[m]]$samples[(initial+1):niter.test[i],])', sep='')))
    }
    eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

    gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
  }
} else {print('error - n_node000 < 1 makes no sense')}

#plot(niter.test, gr.test)

# save
if(DO_WRITE_RDATA) {save.image(file=paste(output_dir,'GEOCARB_MCMC_',appen,'_',today,appen2,'.RData', sep=''))}
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

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

if(FALSE) {#==========================

acf_cutoff <- 0.05
lag_max <- 0.01*niter_mcmc # if we end up with fewer than 100 samples, what are we even doing?
niter_thin <- rep(0, nmodel); names(niter_thin) <- types.of.model
for (model in types.of.model) {
  for (p in 1:length(parnames_all[[model]])) {
    if(n_node000 > 1) {acf_tmp <- acf(x=chains_burned[[model]][[1]][,p], plot=FALSE, lag.max=lag_max)}
    else {acf_tmp <- acf(x=chains_burned[[model]][,p], plot=FALSE, lag.max=lag_max)}
    niter_thin[[model]] <- max(niter_thin[[model]], acf_tmp$lag[which(acf_tmp$acf < acf_cutoff)[1]])
  }
  nthin <- max(niter_thin, na.rm=TRUE)
  if(n_node000 > 1) {
    for (m in 1:n_node000) {
      chains_burned_thinned[[model]][[m]] <- chains_burned[[model]][[m]][seq(from=1, to=nrow(chains_burned[[model]][[m]]), by=nthin),]
    }
  } else {
    chains_burned_thinned[[model]] <- chains_burned[[model]][seq(from=1, to=nrow(chains_burned[[model]]), by=nthin),]
  }
}

}#====================================


# thin to a target number of samples?
if(TRUE) {#===========================

n.sample <- 100000

if(n_node000 == 1) {
  ind.sample <- sample(x=1:nrow(chains_burned), size=n.sample, replace=FALSE)
  chains_burned_thinned <- chains_burned[ind.sample,]
} else {
  n.sample.sub <- rep(NA, n_node000)
  # the case where desired sample size is divisible by the number of chains
  if(round(n.sample/n_node000) == n.sample/n_node000) {
    n.sample.sub[1:n_node000] <- n.sample/n_node000
  } else {
  # the case where it is not
    n.sample.sub[2:n_node000] <- round(n.sample/n_node000)
    n.sample.sub[1] <- n.sample - sum(n.sample.sub[2:n_node000])
  }
  for (m in 1:n_node000) {
    ind.sample <- sample(x=1:nrow(chains_burned[[m]]), size=n.sample.sub[m], replace=FALSE)
    chains_burned_thinned[[m]] <- chains_burned[[m]][ind.sample,]
  }
}

}#====================================


# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).

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

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
if(DO_WRITE_NETCDF) {

lmax=0
for (i in 1:length(parnames_calib)){lmax=max(lmax,nchar(parnames_calib[i]))}

## Name the file
filename.parameters = paste(output_dir,'geocarb_calibratedParameters_',appen,'_',today,appen2,'.nc',sep="")

dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('geocarb_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
covjump.var <- ncvar_def('covjump', '', list(dim.parameters,dim.parameters), -999)
outnc <- nc_create(filename.parameters, list(parameters.var,parnames.var, covjump.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames_calib)
ncvar_put(outnc, covjump.var, covjump.posterior)
nc_close(outnc)

}
##==============================================================================


##==============================================================================
## End
##==============================================================================
