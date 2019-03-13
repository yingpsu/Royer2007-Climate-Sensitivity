##==============================================================================
## GEOCARB_sensitivity_driver.R
##
## Sobol' sensitivity experiment with GEOCARBSULFvolc model
##
## Builds off a previous calibration ensemble output for analysis
## (posterior samples)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

## Clear workspace
rm(list=ls())

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 6000000
.Nboot <- 10000
appen <- 'priors'
.confidence <- 0.9 # for bootstrap CI
.second <- FALSE    # calculate second-order indices?
l_parallel <- FALSE # use parallel evaluation of ensembles in Sobol' integration?
do_sample_tvq <- TRUE
do_precal <- TRUE
USE_LENTON_FSR <- FALSE
USE_ROYER_FSR <- TRUE

# latin hypercube precalibration
sens='NS' # valid values:  L1, L2, NS, pres

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_unc.csv'
#filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_tvq_all-const.csv'
#filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all-const.csv'
filename.calibout <- '../output/geocarb_calibratedParameters_unc_26Feb2019sn.nc'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB/R')
  .Ncore <- 2
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 15  # use multiple cores to process large data?
}

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

##==============================================================================
## Data
##=====

source('GEOCARB-2014_getData.R')

co2_uncertainty_cutoff <- 20

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

# redo with reduced data set

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
if(do_sample_tvq) {
  source('GEOCARB-2014_parameterSetup_tvq.R')
  source('model_forMCMC_tvq.R')
} else {
  source('GEOCARB-2014_parameterSetup.R')
  source('model_forMCMC.R')
}
n_parameters <- length(parnames_calib)
##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distribution bounds
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
# the sensitivity analysis call below.
##==============================================================================


##==============================================================================
## Set up
##=======

# need the physical model
#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `var` uncertainty statistical parameter

# needed packages
#install.packages('sensitivity')
#install.packages('sn')
#install.packages('foreach')
#install.packages('doParallel')
#install.packages('ncdf4')
library(sensitivity)
library(sn)
library(foreach)
library(doParallel)
library(ncdf4)
library(lhs)

##==============================================================================



##==============================================================================

## Read calibration results file, separate into two parameter data frames
ncdata <- nc_open(filename.calibout)
parameters_mcmc <- t(ncvar_get(ncdata, 'geocarb_parameters'))
nc_close(ncdata)
n_ensemble <- nrow(parameters_mcmc)


if(do_precal) {


################
### LHS?


n_sample <- 2*n_sample
n_parameters <- length(parnames_calib)

## draw parameters by Latin Hypercube (Sample)
parameters_lhs <- randomLHS(n_sample, n_parameters)

## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
colnames(par_calib) <- parnames_calib
for (i in 1:n_const_calib) {
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    par_calib[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    par_calib[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}
for (i in (n_const_calib+1):length(parnames_calib)) {
  par_calib[,i] <- qbeta(p=parameters_lhs[,i], shape1=5, shape2=5)
}

tbeg <- proc.time()

model_out <- sapply(1:n_sample, function(ss) {
                    model_forMCMC(par_calib=par_calib[ss,],
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
                                  do_sample_tvq=do_sample_tvq,
                                  par_time_center=par_time_center,
                                  par_time_stdev=par_time_stdev)[,'co2']})

ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < 0) |
      any(model_out[,ss] > 1e4) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}

parameters_good <- par_calib[-ibad,]
colnames(parameters_good) <- parnames_calib

tend <- proc.time()

# report success rate
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))
print(paste('with success rate of ',nrow(parameters_good),'/',n_sample, sep=''))

# save precalibration ensemble
today=Sys.Date(); today=format(today,format="%d%b%Y")
saveRDS(parameters_good, paste('../output/precal_parameters_',today,'.rds', sep=''))


indAvailable <- 1:nrow(parameters_good)
indA <- sample(indAvailable, size=floor(length(indAvailable)/2), replace=FALSE)
indAvailable <- indAvailable[-indA]
indB <- sample(indAvailable, size=length(indA), replace=FALSE)
parameters_sampleA <- parameters_good[indA,]
parameters_sampleB <- parameters_good[indB,]
colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib

################


} else {


## Sample parameters (need 2 data frames of size n_sample each)
if(n_ensemble < 2*n_sample) {print('ERROR - not enough simulations for 2 samples of requested size')}
indAvailable <- 1:n_ensemble
indA <- sample(indAvailable, size=n_sample, replace=FALSE)
indAvailable <- indAvailable[-indA]
indB <- sample(indAvailable, size=n_sample, replace=FALSE)
parameters_sampleA <- parameters_mcmc[indA,]
parameters_sampleB <- parameters_mcmc[indB,]
colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib
##if (ncol(parameters_mcmc)==68) {colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib
##} else if(ncol(parameters_mcmc)==56) {colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_const}


}


##==============================================================================

## Get a reference simulation for integrated sensitivity measure (if using L1, e.g.)
model_out <- model_forMCMC(par_calib=par_calib0,
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
                           do_sample_tvq=do_sample_tvq,
                           par_time_center=par_time_center,
                           par_time_stdev=par_time_stdev)
model_ref <- model_out[,'co2']
model_age <- model_out[,'age']

l_scaled <- TRUE
export_names <- c('model_forMCMC', 'run_geocarbF', 'age', 'ageN',
                  'par_fixed0', 'parnames_calib', 'parnames_fixed', 'parnames_time',
                  'ind_const_calib', 'ind_time_calib', 'ind_const_fixed',
                  'ind_time_fixed', 'input', 'ind_expected_time',
                  'ind_expected_const', 'iteration_threshold', 'l_scaled', 'sens',
                  'model_ref', 'data_calib', 'ind_mod2obs',
                  'do_sample_tvq', 'par_time_center', 'par_time_stdev')

n_boot <- .Nboot
conf <- .confidence
Ncore <- .Ncore

## Actually run the Sobol'

source('sobolTony.R')

if (!l_parallel) {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, parnames_time=parnames_time,
                     age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const, model_ref=model_ref,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     do_sample_tvq=do_sample_tvq, par_time_center=par_time_center, par_time_stdev=par_time_stdev,
                     n_boot=n_boot, conf=conf, second=.second))
} else {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, parnames_time=parnames_time,
                     age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const, model_ref=model_ref,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     do_sample_tvq=do_sample_tvq, par_time_center=par_time_center, par_time_stdev=par_time_stdev,
                     parallel=TRUE, n_core=Ncore, export_names=export_names,
                     n_boot=n_boot, conf=conf, second=.second))
}

print(paste('Sobol simulations took ',t.out[3],' seconds', sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.sobol <- paste('../output/sobol_sens',sens,'_',appen,'_',today,'.rds', sep='')
saveRDS(s.out, filename.sobol)

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T)
max_conf_int <- max(max(s.out$S[,3]-s.out$S[,2]), max(s.out$T[,3]-s.out$T[,2]))
print(paste('max. conf int=',max_conf_int, ' want less than:  0.1 * max. sensitivity index=',max_sens_ind, sep=''))
##==============================================================================



##==============================================================================
## Write indices, results we'd need to file

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_sens',sens,'_',appen,'_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_sens',sens,'_',appen,'_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames_calib, s.out$S, s.out$T))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

if (.second) {
  headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                              'S2_conf_high'), nrow=1)

  output.2nd <- data.frame(cbind( s.out$S2.names,s.out$S2 ))
  write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
  write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
}
##==============================================================================



##==============================================================================
## End
##==============================================================================
