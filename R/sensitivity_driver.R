##==============================================================================
## sensitivity_driver.R
##
## Sobol' sensitivity experiment with GEOCARBSULFvolc model
##
## Builds off a previous calibration ensemble output for analysis
## (posterior samples)
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

## Clear workspace
rm(list=ls())

## Set testing number of samples and file name appendix here
## if there aren't enough samples on the MCMC output file, will break.
n_sample <- 1000000 # note this will be doubled if using LHS precalibration
.Nboot <- 1000
appen <- 'test'
.confidence <- 0.9 # for bootstrap CI
.second <- TRUE    # calculate second-order indices?

DO_PARALLEL <- TRUE # use parallel evaluation of ensembles in Sobol' integration?
DO_SAMPLE_TVQ <- TRUE  # sample time series uncertainty by CDF parameters?
USE_LENTON_FSR <- FALSE
USE_ROYER_FSR <- TRUE

# latin hypercube precalibration
sens='NS' # valid values:  L1, L2, NS, pres

# calibration parameters and input data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_unc.csv'

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB/R')
  .Ncore <- 2
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('~/work/codes/GEOCARB/R')
  .Ncore <- 6  # use multiple cores to process large data?
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
##==============================================================================


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
source('parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')
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
#source('run_geocarbF.R')n   # version without the extra uncertainty parameter
source('run_geocarbF_unc.R') # version with extra `stdev` uncertainty statistical parameter

# needed packages
library(sensitivity)
library(sn)
library(foreach)
library(doParallel)
library(lhs)

##==============================================================================



##==============================================================================
## Draw parameter samples
##=======================

## double the sample size, to generate two independent arrays (Sobol, 2001)
n_sample <- 2*n_sample

## scale up to the actual parameter distributions

## draw parameters by Latin Hypercube (Sample)
set.seed(2019)
parameters_lhs <- randomLHS(n_sample, n_parameters)
par_calib <- parameters_lhs  # initialize

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
##==============================================================================



##==============================================================================
## parameter precalibration
##=========================

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
                                  do_sample_tvq=DO_SAMPLE_TVQ,
                                  par_time_center=par_time_center,
                                  par_time_stdev=par_time_stdev)[,'co2']})
ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < .lower_bound_co2) |
      any(model_out[,ss] > .upper_bound_co2) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}
parameters_good <- par_calib[-ibad,]
tend <- proc.time()
# report success rate
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))
print(paste('with success rate of ',nrow(parameters_good),'/',n_sample, sep=''))
# save precalibration ensemble
today=Sys.Date(); today=format(today,format="%d%b%Y")
saveRDS(parameters_good, paste('../output/precal_parameters_',today,'.rds', sep=''))

set.seed(9102)
colnames(parameters_good) <- parnames_calib
indAvailable <- 1:nrow(parameters_good)
indA <- sample(indAvailable, size=floor(length(indAvailable)/2), replace=FALSE)
indAvailable <- indAvailable[-indA]
indB <- sample(indAvailable, size=length(indA), replace=FALSE)
parameters_sampleA <- parameters_good[indA,]
parameters_sampleB <- parameters_good[indB,]
colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib
##==============================================================================



##==============================================================================
## Sobol analysis
##===============

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
                           do_sample_tvq=DO_SAMPLE_TVQ,
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
                  'DO_SAMPLE_TVQ', 'par_time_center', 'par_time_stdev')

n_boot <- .Nboot
conf <- .confidence
Ncore <- .Ncore

## Actually run the Sobol'

source('sobolTony.R')

if (!DO_PARALLEL) {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, parnames_time=parnames_time,
                     age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const, model_ref=model_ref,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center,
                     par_time_stdev=par_time_stdev, n_boot=n_boot, conf=conf, second=.second))
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
                     do_sample_tvq=DO_SAMPLE_TVQ, par_time_center=par_time_center,
                     par_time_stdev=par_time_stdev, parallel=TRUE, n_core=Ncore,
                     export_names=export_names, n_boot=n_boot, conf=conf, second=.second))
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
  headers.2nd <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
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
