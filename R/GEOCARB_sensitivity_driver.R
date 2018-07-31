##==============================================================================
## precal_sens_sobolOnly.R
##
## Precalibration and
## sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Builds off a previous LHS precalibration sample
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

## Clear workspace
rm(list=ls())

## Set testing number of samples and file name appendix here
n_sample <- 100
appen <- 'TEST'
.Nboot <- 100
.confidence <- 0.9 # for bootstrap CI
.scheme <- 'A' # A = first and total indices; B = first, second and total
l_parallel <- TRUE

co2_uncertainty_cutoff <- 20

# latin hypercube precalibration
alpha <- 0
sens='NS'

filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB/R')
  .Ncore <- 2
  filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL1_30Mar2018.csv'
  #filename_in <- '../output/geocarb_precalibration_parameters_alpha10_sensL2_25Mar2018.csv'
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 15  # use multiple cores to process large data?
  filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL1_01Apr2018.csv'
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

filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_GAMMA-co2_31Jul2018.csv'
#filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_LN-co2_31Jul2018.csv'
#filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_06Jun2017.csv'
source('GEOCARB-2014_getData.R')

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
  if(length(ind_remove) > 0) {data_calib <- data_calib[-ind_remove,]}
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
n_parameters <- length(parnames_calib)
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
bounds_calib <- mat.or.vec(nr=n_parameters, nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:n_parameters) {
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

# needed packages
#install.packages('sensitivity')
#install.packages('sn')
#install.packages('lhs')
#install.packages('foreach')
#install.packages('doParallel')
library(sensitivity)
library(sn)
library(lhs)
library(foreach)
library(doParallel)
##==============================================================================


## this is where the LHS precalibration would go...


##==============================================================================

## Read KDE results file, separate into parameters and the bandwidths
## (Moved to beginning becase this depends on remote vs local, whehter small testing
##  or ready to go large parameter sets)
#filename_in <- filename_out
#alpha <- 0; filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL1_30Mar2018.csv'
#alpha <- 0; filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL1_01Apr2018.csv'
#alpha <- 0; filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL2_24Mar2018.csv'
#alpha <- 0.10; filename_in <- '../output/geocarb_precalibration_parameters_alpha10_sensL2_25Mar2018.csv'
#alpha <- 0.34; filename_in <- '../output/geocarb_precalibration_parameters_alpha34_sensL2_24Mar2018.csv'
parameters_node <- read.csv(filename_in)
n_node <- nrow(parameters_node)-1
bandwidths <- parameters_node[n_node+1,]
parameters_node <- parameters_node[-(n_node+1),]

##==============================================================================

## Get a reference simulation for integrated sensitivity measure (if using L1, e.g.)
model_ref <- model_forMCMC(par_calib=par_calib0,
              par_fixed=par_fixed0,
              parnames_calib=parnames_calib,
              parnames_fixed=parnames_fixed,
              age=age,
              ageN=ageN,
              ind_const_calib=ind_const_calib,
              ind_time_calib=ind_time_calib,
              ind_const_fixed=ind_const_fixed,
              ind_time_fixed=ind_time_fixed,
              ind_expected_time=ind_expected_time,
              ind_expected_const=ind_expected_const,
              iteration_threshold=iteration_threshold)[,'co2']

l_scaled <- TRUE
export_names <- c('model_forMCMC', 'run_geocarbF',
                  'par_fixed0', 'parnames_calib', 'parnames_fixed', 'age', 'ageN',
                  'ind_const_calib', 'ind_time_calib', 'ind_const_fixed',
                  'ind_time_fixed', 'input', 'ind_expected_time',
                  'ind_expected_const', 'iteration_threshold', 'l_scaled', 'sens',
                  'model_ref', 'data_calib', 'ind_mod2obs')

n_boot <- .Nboot
conf <- .confidence
Ncore <- .Ncore

## Sample parameters (need 2 data frames) by taking directly from the precalibration?
n_half <- floor(0.5*nrow(parameters_node))

## If n_sample is given, then use parameter samples of that size
if (n_sample > 0) {n_half <- n_sample}
parameters_sampleA <- parameters_node[1:n_half,]
parameters_sampleB <- parameters_node[(n_half+1):(2*n_half),]
colnames(parameters_sampleA) <- colnames(parameters_sampleB) <- parnames_calib

## Actually run the Sobol'

source('sobolTony.R')

if (!l_parallel) {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time,
                     ind_expected_const=ind_expected_const,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     n_boot=n_boot, conf=conf))
} else {
  t.out <- system.time(s.out <- sobolTony(parameters_sampleA, parameters_sampleB, sens,
                     par_fixed=par_fixed0, parnames_calib=parnames_calib,
                     parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                     ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                     ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                     input=input, ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                     iteration_threshold=iteration_threshold, data_calib=data_calib,
                     parallel=TRUE, n_core=Ncore, export_names=export_names,
                     n_boot=n_boot, conf=conf))
}

print(paste('Sobol simulations took ',t.out[3],' seconds', sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.sobol <- paste('../output/sobol_alpha',100*alpha,'_',appen,'_',today,'.rds', sep='')
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
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_alpha',100*alpha,'_sens',sens,'_',appen,'_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_alpha',100*alpha,'_sens',sens,'_',appen,'_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames_calib, s.out$S, s.out$T))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

if (.scheme=='B') {
  headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                              'S2_conf_high'), nrow=1)
  output2.indices <- s.out$S2[,1]
  output2.conf1   <- s.out$S2[,4]
  output2.conf2   <- s.out$S2[,5]

  # 2nd order index names ordered as: (assuming 39 parameters)
  # 1. parnames.sobol[1]-parnames.sobol[2]
  # 2. parnames.sobol[1]-parnames.sobol[3]
  # 3. parnames.sobol[1]-parnames.sobol[4]
  # ... etc ...
  # 38. parnames.sobol[1]-parnames.sobol[39] << N=2:39 => p1-p[N]
  # 39. parnames.sobol[2]-parnames.sobol[3]
  # 40. parnames.sobol[2]-parnames.sobol[4]
  # 38+37. parnames.sobol[2]-parnames.sobol[39] << N=3:39 => p2-p[N]
  # ... etc ...
  names2  <- rownames(s.out$S2)
  names2a <- rep(NA, length(names2))
  names2b <- rep(NA, length(names2))
  cnt <- 1
  for (i in seq(from=1, to=(length(parnames_calib)-1), by=1)) {         # i = index of first name
      for (j in seq(from=(i+1), to=(length(parnames_calib)), by=1)) {   # j = index of second name
          names2a[cnt] <- parnames_calib[i]
          names2b[cnt] <- parnames_calib[j]
          cnt <- cnt+1
      }
  }

  output.2nd <- data.frame(cbind( names2a,
                                  names2b,
                                  output2.indices,
                                  output2.conf1,
                                  output2.conf2 ))
  write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
  write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
              quote=FALSE    , row.names = FALSE , col.names=FALSE)
}
##==============================================================================

if(FALSE){



}


##==============================================================================
## End
##==============================================================================
