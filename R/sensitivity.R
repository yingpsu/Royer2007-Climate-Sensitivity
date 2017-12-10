##==============================================================================
## sensitivity_geocarb.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())

co2_uncertainty_cutoff <- 20


## which experiments to do?
ldosobol <- TRUE
ldomorris <- FALSE
ldoparallel <- TRUE


if(Sys.info()['nodename']=='Tonys-MacBook-Pro.local') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
  .Ncore <- 0
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 5  # use multiple cores to process large data?
}



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
if(co2_uncertainty_cutoff > 0) {
  co2_halfwidth <- 0.5*(data_calib$co2_high - data_calib$co2_low)
  ind_filter <- which(co2_halfwidth < co2_uncertainty_cutoff)
  data_calib <- data_calib[-ind_filter,]
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
##==============================================================================


##==============================================================================
## Run the sensitivity analysis
##=============================

#install.packages('sensitivity')
#install.packages('sn')
library(sensitivity)
library(sn)

## Method of Morris
## Send in bounds of [0,1] for each parameter, and then scale within the
## `log_like_sensitivity` function to the parameters' distributions (using CDF)

repetitions <- 20
levels <- 500
grid_jump <- round(levels/2)
perc_lower <- rep(0.01, n_parameters)
perc_upper <- rep(0.99, n_parameters)

if(ldomorris) {
tbeg <- proc.time()
mom <- morris(model=log_like_sensitivity, factors=parnames_calib, r=repetitions,
              binf=perc_lower, bsup=perc_upper, scale=FALSE,
              design=list(type="oat", levels=levels, grid.jump=grid_jump),
              par_fixed=par_fixed0, parnames_calib=parnames_calib,
              parnames_fixed=parnames_fixed, age=age, ageN=ageN,
              ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
              ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
              input=input,
              data_calib=data_calib, ind_mod2obs=ind_mod2obs,
              ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
              iteration_threshold=iteration_threshold)
tend <- proc.time()
print(paste(length(mom$y),' simulations took ',(tend-tbeg)[3]/60,' minutes',sep=''))

mu <- apply(mom$ee, 2, mean)
mu_star <- apply(mom$ee, 2, function(x) mean(abs(x)))
sigma <- apply(mom$ee, 2, sd)
stderr_mu <- sigma/sqrt(repetitions)
mom_table <- cbind(mu, mu_star, sigma)
}
##==============================================================================


##==============================================================================
## Sobol' Method

#install.packages('lhs')
library(lhs)

## Sobol' wrapper - assumes uniform distributions on parameters
geocarb_sobol_ser <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, data_calib, ind_mod2obs, ind_expected_time,
                          ind_expected_const, iteration_threshold) {
  ll_output <- log_like_sensitivity(par_calib_scaled,
                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input,
                data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold)
  ll_output_centered <- ll_output - mean(ll_output)
  return(ll_output_centered)
}

## Sobol' wrapper - assumes uniform distributions on parameters
geocarb_sobol_par <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, data_calib, ind_mod2obs, ind_expected_time,
                          ind_expected_const, iteration_threshold, Ncore) {

  #install.packages('foreach')
  #install.packages('doParallel')
  library(foreach)
  library(doParallel)
  cores=detectCores()
#  cl <- makeCluster(cores[1]-1) #not to overload your computer
  cl <- makeCluster(Ncore)
  print(paste('Starting cluster with ',Ncore,' cores', sep=''))
  registerDoParallel(cl)

  ## initialize output/input
  n_simulations <- nrow(par_calib_scaled)
  export_names <- c('log_like_sensitivity', 'model_forMCMC', 'run_geocarbF',
                    'par_fixed0','parnames_calib','parnames_fixed','age','ageN','ind_const_calib','ind_time_calib','ind_const_fixed','ind_time_fixed','input','data_calib','ind_mod2obs','ind_expected_time','ind_expected_const','iteration_threshold')
  output <- rep(0, n_simulations)

  finalOutput <- foreach(ii=1:n_simulations, .combine=c,
                                   .packages=c('sn','stats'),
                                   .export=export_names,
                                   .inorder=FALSE) %dopar% {

    dyn.load("../fortran/run_geocarb.so")
    ll_output <- log_like_sensitivity(as.numeric(par_calib_scaled[ii,]),
                    par_fixed=par_fixed0, parnames_calib=parnames_calib,
                    parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                    ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                    ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                    input=input,
                    data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                    ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                    iteration_threshold=iteration_threshold)
    output[ii] <- ll_output
  }
  stopCluster(cl)
  ll_output_centered <- output - mean(output)
  return(ll_output_centered)
}

# getting about 290000 simulations in 2 minutes (# simulations = (p+2)*n_sample)
n_sample <- 200
n_bootstrap <- 100

## Sample parameters (need 2 data frames)
parameters_lhs1 <- randomLHS(n_sample, n_parameters)
parameters_lhs2 <- randomLHS(n_sample, n_parameters)

## Trim so you aren't sampling the extreme cases?
alpha <- 0.5
parameters_lhs1 <- (1-alpha)*parameters_lhs1 + 0.5*alpha
parameters_lhs2 <- (1-alpha)*parameters_lhs2 + 0.5*alpha

## Need data frames as input
parameters_lhs1 <- data.frame(parameters_lhs1)
parameters_lhs2 <- data.frame(parameters_lhs2)
colnames(parameters_lhs1) <- colnames(parameters_lhs2) <- parnames_calib

## Actually run the Sobol'
if(ldosobol) {
  if(ldoparallel ) {
    t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_par,
                             parameters_lhs1,
                             parameters_lhs2,
                             scheme='A',
                             nboot=n_bootstrap,
                             par_fixed=par_fixed0, parnames_calib=parnames_calib,
                             parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                             ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                             ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                             input=input,
                             data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                             ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                             iteration_threshold=iteration_threshold, Ncore=.Ncore))
  } else {
    t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_ser,
                             parameters_lhs1,
                             parameters_lhs2,
                             scheme='A',
                             nboot=n_bootstrap,
                             par_fixed=par_fixed0, parnames_calib=parnames_calib,
                             parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                             ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                             ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                             input=input,
                             data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                             ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                             iteration_threshold=iteration_threshold))
  }
## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('max. sensitivity index=',max_sens_ind,' // max. conf int=',max_conf_int,sep=''))
}

today <- Sys.Date(); today <- format(today,format="%d%b%Y")
save.image(file=paste('geocarb_sensitivity_',today,'.RData', sep=''))

##==============================================================================



##==============================================================================
## End
##==============================================================================
