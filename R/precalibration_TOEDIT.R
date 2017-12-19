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

n_sample <- 2e5


if(Sys.info()['nodename']=='Tonys-MacBook-Pro.local') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
  .Ncore <- 2
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 15  # use multiple cores to process large data?
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


##==============================================================================
## Draw parameters by latin hypercube
##===================================

parameters_lhs <- randomLHS(n_sample, n_parameters)


## Trim so you aren't sampling the extreme cases?
alpha <- .5
parameters_lhs <- (1-alpha)*parameters_lhs + 0.5*alpha


## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
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
##==============================================================================


##==============================================================================
## Run precalibration simulations
##===============================

model_out <- sapply(1:n_sample, function(ss) {
    model_forMCMC(par_calib=par_calib[ss,],
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
                  iteration_threshold=iteration_threshold)[,'co2']})

#source('precalibration_doruns.R')

ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < 0) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}

parameters_good <- par_calib[-ibad,]
colnames(parameters_good) <- parnames_calib
##==============================================================================


##==============================================================================
## Fit KDEs to sample from for each parameter
##===========================================

pdf_all <- vector('list', n_parameters)
for (pp in 1:n_parameters){
  tmp <- density(parameters_good[,pp], kernel='gaussian', n=100)
  pdf_all[[pp]] <- tmp; names(pdf_all)[pp] <- parnames_calib[pp]
}

## Write a CSV file with the successful parameter combinations and bandwidths
bandwidths <- rep(NA,n_parameters)
for (pp in 1:n_parameters){
  bandwidths[pp] <- pdf_all[[pp]]$bw
}
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_out <- paste('geocarb_precalibration_parameters_',today,'.csv', sep='')
write.csv(rbind(parameters_good, bandwidths), file=filename_out, row.names=FALSE)
##==============================================================================


##==============================================================================
## Function for sampling from KDEs (inflating a LHS?)
##===================================================

## Read KDE results file, separate into parameters and the bandwidths
filename_in <- 'geocarb_precalibration_parameters_19Dec2017.csv'
parameters_node <- read.csv(filename_in)
n_node <- nrow(parameters_node)-1
bandwidths <- parameters_node[n_node+1,]
parameters_node <- parameters_node[-(n_node+1),]

kde_sample <- function(n_sample, nodes, bandwidths) {
  # preliminaries
  n_node <- nrow(nodes)
  n_par <- length(bandwidths)
  if(n_sample > n_node) {print('ERROR: n_sample > n_node')}

  # choose the node rows out of array `nodes`
  i_sample <- sample(x=1:n_node, size=n_sample, replace=FALSE)

  # sample normal random from around each of the parameter nodes from that row
  # (this achieves joint sampling)
  par_sample <- t(sapply(1:n_sample, function(i) {
       rnorm(n=n_parameters, mean=as.numeric(parameters_node[i,]), sd=as.numeric(bandwidths))}))

  return(par_sample)
}
##==============================================================================



##==============================================================================

## Sobol' wrapper - assumes uniform distributions on parameters
geocarb_sobol_co2_ser <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, ind_expected_time,
                          ind_expected_const, iteration_threshold) {
  co2_output <- sensitivity_co2(par_calib_scaled,
                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input, ind_expected_time=ind_expected_time,
                ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold)
#  co2_output_centered <- co2_output - mean(co2_output, na.rm=TRUE)
  co2_output_centered <- co2_output - mean(co2_output[is.finite(co2_output)])
  co2_output_centered[is.infinite(co2_output)] <- 0
  return(co2_output_centered)
}

n_sample <- 100
n_bootstrap <- 5
Ncore <- 1

source('GEOCARB_sensitivity_co2.R')

## Sample parameters (need 2 data frames)
parameters_sample1 <- kde_sample(n_sample, parameters_node, bandwidths)
parameters_sample2 <- kde_sample(n_sample, parameters_node, bandwidths)
colnames(parameters_sample1) <- colnames(parameters_sample2) <- parnames_calib

## Actually run the Sobol'
t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_co2_ser,
                           parameters_sample1,
                           parameters_sample2,
                           scheme='A',
                           nboot=n_bootstrap,
                           par_fixed=par_fixed0, parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                           ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                           input=input,
                           ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold))

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('0.1 * max. sensitivity index=',max_sens_ind,' // max. conf int=',max_conf_int,sep=''))

##==============================================================================



##==============================================================================
## Try Morris?

alpha <- 0.1
perc_lower <- rep(0.5*alpha  , n_parameters)
perc_upper <- rep(1-0.5*alpha, n_parameters)

repetitions <- c(20, 40, 80, 160, 320, 640, 1280)
levels <- c(20, 40, 80, 160, 320, 640, 1280)

n_rep <- length(repetitions)
n_lev <- length(levels)

names_rep <- paste('r',repetitions[1],sep='')
if(n_rep>1) {for (rr in 2:n_rep) {names_rep <- c(names_rep, paste('r',repetitions[rr], sep=''))}}
names_lev <- paste('l',levels[1],sep='')
if(n_lev>1) {for (ll in 2:n_lev) {names_lev <- c(names_lev, paste('l',repetitions[ll], sep=''))}}

mom <- vector('list', n_rep); names(mom) <- names_rep
for (rr in names_rep) {
  mom[[rr]] <- vector('list', n_lev); names(mom[[rr]]) <- names_lev
}

## initialize other output
mu <- mu_star <- med <- med_star <- sigma <- stderr_mu <- mom_table <- names_signif <- mom

for (rr in 1:n_rep) {
  for (ll in 1:n_lev) {
    tbeg <- proc.time()
    grid_jump <- round(0.5*levels[ll])
    mom[[rr]][[ll]] <- morris(model=sensitivity_co2, factors=parnames_calib, r=repetitions[rr],
                              binf=perc_lower, bsup=perc_upper, scale=FALSE,
                              design=list(type="oat", levels=levels[ll], grid.jump=grid_jump),
                              par_fixed=par_fixed0, parnames_calib=parnames_calib,
                              parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                              ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                              ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                              input=input,
                              ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                              iteration_threshold=iteration_threshold, l_scaled=FALSE)
    tend <- proc.time()
    print(paste(length(mom[[rr]][[ll]]$y),' simulations took ',(tend-tbeg)[3]/60,' minutes',sep=''))
    mu[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, mean)
    mu_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) mean(abs(x)))
    med[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, median)
    med_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) median(abs(x)))
    sigma[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, sd)
    stderr_mu[[rr]][[ll]] <- sigma[[rr]][[ll]]/sqrt(repetitions[rr])
    mom_table[[rr]][[ll]] <- cbind(mu[[rr]][[ll]], mu_star[[rr]][[ll]], sigma[[rr]][[ll]], stderr_mu[[rr]][[ll]])
  }
}
## Get the set of significant (mu_star > 2*stderr_mu) parameters for each of the
## simulations; take the intersection (union?) of the "okay" ones as the
## parameter set for calibration?
names_all <- NULL
for (rr in 1:n_rep) {
  for (ll in 1:n_lev) {
    names_signif[[rr]][[ll]] <- parnames_calib[which(mu_star[[rr]][[ll]] > 2*stderr_mu[[rr]][[ll]])]
    names_all <- union(names_all, names_signif[[rr]][[ll]])
  }
}

##==============================================================================



##==============================================================================


##==============================================================================
## Save progress
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
save.image(file=paste('geocarb_precalibration_',today,'.RData', sep=''))


##==============================================================================



##==============================================================================
## End
##==============================================================================
