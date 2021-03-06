##==============================================================================
## mwg_driver.R
##
## Driver for the Metropolis-within-Gibbs sampling.
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


rm(list=ls())


library(compiler)
enableJIT(3)

niter_mcmc000 <- 5e4   # number of MCMC iterations per node (Markov chain length)
do_adapt000 <- TRUE    # adapt Metropolis proposal variances?
n_node000 <- 1         # number of CPUs to use
#appen <- 'sig18+GLAC+LIFE'
appen <- 'sig18'
#appen <- 'all'
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
co2_uncertainty_cutoff <- 20

DO_INIT_UPDATE <- TRUE

DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- FALSE

filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.par_fixed  <- '../output/par_deoptim_OPT1_04Jul2018.rds'
filename.covariance <- paste('../output/par_LHS2_',appen,'_04Jul2018.RData', sep='')

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('~/codes/GEOCARB/R')
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('~/codes/GEOCARB/R')
}

library(sn)
library(ncdf4)

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
  data_calib <- data_calib[-ind_filter,]    # removing all of the possibly troublesome points
  ##data_calib <- data_calib[-ind_remove,]    # remove only those the revised range does not help
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
  ##sd <- 2.4*2.4/length(par_calib0)
  sd <- 2.4*2.4/1  # only /1 because metropolis-within-gibbs (MWG) sampling
  eps <- 0.0001
  step_mcmc <- sd*cov(parameters_good) + sd*eps*diag(x=1, length(par_calib0))
}
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

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# need the likelihood function and prior distributions, and the MWG sampler
source('GEOCARB-2014_calib_likelihood.R')
source('mwg.R')

# set up and run the actual calibration
niter_mcmc <- niter_mcmc000
do_adapt <- do_adapt000
stopadapt_mcmc <- round(niter_mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)

##==============================================================================
## Actually run the calibration

tbeg <- proc.time()

if(n_node000==1) {
  mcmc_out <- mcmc_mwg(log_post, n=niter_mcmc, init=par_calib0, scale=diag(step_mcmc),
                     adapt=do_adapt, n.start=max(5000,round(0.2*niter_mcmc)),
                     parallel=FALSE, n.core=1,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                     data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                     ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                     iteration_threshold=iteration_threshold)
  chain1 <- mcmc_out$samples
} else {
  mcmc_out <- mcmc_mwg(log_post, n=niter_mcmc, init=par_calib0, scale=diag(step_mcmc),
                     adapt=do_adapt, n.start=max(5000,round(0.2*niter_mcmc)),
                     parallel=FALSE, n.core=1,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                     data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                     ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                     iteration_threshold=iteration_threshold)
  chain1 <- mcmc_out$samples
}

tend <- proc.time()
print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

if(FALSE) {
par(mfrow=c(2,1))
ics <- match('deltaT2X', parnames_calib)
plot(chain1[,ics], type='l', ylab=parnames_calib[ics], xlab='Iteration')
}

# save
if(DO_WRITE_RDATA) {
  save.image(file=paste(output_dir,'GEOCARB_MCMC-MWG_',appen,'_',today,'.RData', sep=''))
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
