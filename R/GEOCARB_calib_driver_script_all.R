##==============================================================================
## GEOCARB-2014_calib_driver.R
##
## Read CO2 proxy data. Set which data sets you intend to calibrate using.
## Version for running as script on HPC.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')

niter_mcmc000 <- 1e3   # number of MCMC iterations per node (Markov chain length)
n_node000 <- 1         # number of CPUs to use
#appen <- 'sig18+GLAC+LIFE'
#appen <- 'sig18'
appen <- 'tvq_all'
appen2 <- ''
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
co2_uncertainty_cutoff <- 20

# Distribution fit to each data point in the processing step
#dist <- 'ga'
#dist <- 'be'
#dist <- 'ln'
dist <- 'sn'

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

DO_SAMPLE_TVQ <- TRUE  # sample time series uncertainty by CDF parameters?
DO_INIT_UPDATE <- FALSE
DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- FALSE

filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.par_fixed  <- '../output/par_deoptim_OPT1_04Jul2018.rds'
filename.covariance <- paste('../output/par_LHS2_',appen,'_04Jul2018.RData', sep='')

library(adaptMCMC)
library(ncdf4)
library(sn)

##==============================================================================
## Model parameters and setup
##===========================

if(DO_SAMPLE_TVQ) {
  source('GEOCARB-2014_parameterSetup_tvq.R')
} else {
  source('GEOCARB-2014_parameterSetup.R')
}

if (DO_INIT_UPDATE) {
  # initial fixed parameters estimate:
  par_new <- readRDS(filename.par_fixed)
  for (name in names(par_new)) {
    if (name %in% parnames_fixed) {
      par_fixed0[match(name, parnames_fixed)] <- par_new[name]
    }
    if (name %in% parnames_calib) {
      par_calib0[match(name, parnames_calib)] <- par_new[name]
    }
  }
  # need to strip the names or MCMC will break
  names(par_calib0) <- NULL
  names(par_fixed0) <- NULL

  # initial covariance estimate:
  load(filename.covariance)
  sd <- 2.4*2.4/length(par_calib0)
  eps <- 0.0001
  step_mcmc <- sd*cov(parameters_good) + sd*eps*diag(x=1, length(par_calib0))
}
##==============================================================================


##==============================================================================
## Data
##=====

# need the physical model
if(DO_SAMPLE_TVQ) {source('model_forMCMC_tvq.R')
} else            {source('model_forMCMC.R')}
source('run_geocarbF.R')

source('GEOCARB_fit_likelihood_surface.R')
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

# only need to rearrange if DO_SAMPLE_TVQ
if (DO_SAMPLE_TVQ) {
  bound_lower <- bound_lower[match(parnames_calib, as.character(input$parameter))]
  bound_upper <- bound_upper[match(parnames_calib, as.character(input$parameter))]
}

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the calibration MCMC call below.
##==============================================================================


##==============================================================================
## Pad if only one calibration parameter
## (adaptMCMC requires 2 or more)
##======================================
if(length(parnames_calib)==1){
  parnames_calib <- c(parnames_calib, "padding")
  bounds_calib <- rbind(bounds_calib,c(-Inf,Inf))
  rownames(bounds_calib) <- parnames_calib
  par_calib <- c(par_calib, 0)
  step_mcmc <- c(step_mcmc, 1)
}
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
    tbeg <- proc.time()
    amcmc_out1 <- MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                    scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=max(5000,round(0.05*niter_mcmc)),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                    ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                    ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                    input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                    data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                    ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                    iteration_threshold=iteration_threshold,
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                    do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
    tend <- proc.time()
    chain1 = amcmc_out1$samples
  } else if(n_node000 > 1) {
    tbeg <- proc.time()
    amcmc.par1 <- MCMC.parallel(log_post, n=niter_mcmc, init=par_calib0, n.chain=n_node000, n.cpu=n_node000,
          					dyn.libs=c('../fortran/run_geocarb.so'),
                    packages=c('sn'),
          					adapt=TRUE, list=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
          					gamma=gamma_mcmc, n.start=max(5000,round(0.05*niter_mcmc)),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                    ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                    ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                    input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                    data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                    ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                    iteration_threshold=iteration_threshold,
                    loglikelihood_smoothed=loglikelihood_smoothed, likelihood_fit=likelihood_fit, idx_data=idx_data,
                    do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center, par_time_stdev=par_time_stdev)
    tend <- proc.time()
  }
} else {
  if(n_node000==1) {
    tbeg <- proc.time()
    amcmc_out1 <- MCMC(log_post, n=niter_mcmc, init=par_calib0, adapt=TRUE, acc.rate=accept_mcmc,
                    scale=step_mcmc, gamma=gamma_mcmc, list=TRUE, n.start=max(5000,round(0.05*niter_mcmc)),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, age=age, ageN=ageN,
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
    amcmc.par1 <- MCMC.parallel(log_post, n=niter_mcmc, init=par_calib0, n.chain=n_node000, n.cpu=n_node000,
          					dyn.libs=c('../fortran/run_geocarb.so'),
                    packages=c('sn'),
          					adapt=TRUE, list=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
          					gamma=gamma_mcmc, n.start=max(5000,round(0.05*niter_mcmc)),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, age=age, ageN=ageN,
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

if(FALSE) {
par(mfrow=c(2,1))
ics <- match('deltaT2X', parnames_calib)
plot(chain1[,ics], type='l', ylab=parnames_calib[ics], xlab='Iteration')
}

# save
if(DO_WRITE_RDATA) {
  save.image(file=paste(output_dir,'GEOCARB_MCMC-CON_',appen,'_',today,appen2,'.RData', sep=''))
}

## Extend an MCMC chain?
## Extend and run more MCMC samples?
if(FALSE){
niter_extend <- 1e4
tbeg=proc.time()
amcmc_extend1 = MCMC.add.samples(amcmc_out1, niter_extend,
                                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
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

# visual inspection

#TONY TODO
#TONY TODO
# for now only look at climate sensitivity
#ind_cs <- match('deltaT2X',parnames_calib)
#plot(chain1[,ind_cs], type='l')
#par(mfrow=c(7,8))
#for (p in 1:length(parnames_calib)) {plot(chain1[,p], type='l', ylab=parnames_calib[p])}


#TODO

## Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
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
      eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc.par1[[m]]$samples[1:niter.test[i],])', sep='')))
    }
    eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

    gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
  }
} else {print('error - n_node000 < 1 makes no sense')}

#plot(niter.test, gr.test)

# Heidelberger and Welch
# visual inspection

# which parameters is climate sensitivity?
if(FALSE) {
ics <- which(parnames_calib=='deltaT2X')
plot(chain1[,ics], type='l')
}

#parameters.posterior <- chain1
#covjump <- amcmc_out1$cov.jump

##==============================================================================


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
filename.parameters = paste('geocarb_calibratedParameters_',appen,'_',today,'.nc',sep="")

dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('geocarb_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
covjump.var <- ncvar_def('covjump', '', list(dim.parameters,dim.parameters), -999)
outnc <- nc_create(filename.parameters, list(parameters.var,parnames.var, covjump.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames_calib)
ncvar_put(outnc, covjump.var, covjump)
nc_close(outnc)

}
##==============================================================================


##==============================================================================
## End
##==============================================================================
