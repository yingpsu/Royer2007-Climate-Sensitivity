##==============================================================================
## precal_sens.R
##
## Precalibration and
## sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())

co2_uncertainty_cutoff <- 20

# latin hypercube precalibration
alpha <- 0
n_sample <- 5.12e5
n_sample_max <- 2.56e5 # if asking n_sample > n_sample_max, break into subsamples
sens='L1'

filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'

if(Sys.info()['user']=='tony') {
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
  data_calib <- data_calib[-ind_remove,]
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


tbeg <- proc.time()

if (n_sample <= n_sample_max) {

  ##====================================
  ## Draw parameters by latin hypercube
  ##   then
  ## Run precalibration simulations
  ##====================================

  parameters_lhs <- randomLHS(n_sample, n_parameters)

  ## Trim so you aren't sampling the extreme cases?
  #alpha <- 1-.66
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

  ibad <- NULL
  for (ss in 1:n_sample) {
    if( any(model_out[,ss] < 100) |
        any(model_out[,ss] > 1e4) |
        model_out[58,ss] < 280 | model_out[58,ss] > 400) {
      ibad <- c(ibad,ss)
    }
  }

} else {

  ## Break up into smaller LHS if asking for more than n_sample_max samples

  n_full <- floor(n_sample/n_sample_max)            # full n_sample_max samples
  n_part <- ceiling(n_sample/n_sample_max) - n_full # partial samples (< n_sample_max)

good_parameters <- vector('list', n_full+n_part)

  print('Breaking up into:')
  print(paste('  ',n_full,' samples of size ',n_sample_max, sep=''))
  print(paste('  ',n_part,' samples of size ',n_sample-n_full*n_sample_max, sep=''))

  # do all the full n_sample_max samples
  par_calib <- mat.or.vec(n_sample, n_parameters)

  for (k in 1:n_full) {

    ## start at NULL and add for each subsample the bad runs
    ibad <- NULL

    parameters_lhs <- randomLHS(n_sample_max, n_parameters)

    ## Trim so you aren't sampling the extreme cases?
    #alpha <- 1-.66
    parameters_lhs <- (1-alpha)*parameters_lhs + 0.5*alpha

    ## scale up to the actual parameter distributions
    n_const_calib <- length(ind_const_calib)
    par_calib_subsample <- parameters_lhs  # initialize
    for (i in 1:n_const_calib) {
      row_num <- match(parnames_calib[i],input$parameter)
      if(input[row_num, 'distribution_type']=='gaussian') {
        par_calib_subsample[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      } else if(input[row_num, 'distribution_type']=='lognormal') {
        par_calib_subsample[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
      } else {
        print('ERROR - unknown prior distribution type')
      }
    }
    ind_subsample <- ((k-1)*n_sample_max+1):(k*n_sample_max)
    par_calib[ind_subsample,] <- par_calib_subsample

    # run model, precalibration windows
    model_out <- sapply(1:n_sample_max, function(ss) {
        model_forMCMC(par_calib=par_calib_subsample[ss,],
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
    for (ss in 1:n_sample_max) {
      if( any(model_out[,ss] < 100) |
          any(model_out[,ss] > 1e4) |
          model_out[58,ss] < 280 | model_out[58,ss] > 400) {
        ibad <- c(ibad, ss)
      }
    }
    good_parameters[[k]] <- par_calib_subsample[-ibad,]
  }

  # and a partial sample, if there is one
  if (n_part) {

    ## start at NULL and add for each subsample the bad runs
    ibad <- NULL

    ind_subsample <- (n_full*n_sample_max+1):n_sample

    parameters_lhs <- randomLHS(length(ind_subsample), n_parameters)

    ## Trim so you aren't sampling the extreme cases?
    #alpha <- 1-.66
    parameters_lhs <- (1-alpha)*parameters_lhs + 0.5*alpha

    ## scale up to the actual parameter distributions
    n_const_calib <- length(ind_const_calib)
    par_calib_subsample <- parameters_lhs  # initialize
    for (i in 1:n_const_calib) {
      row_num <- match(parnames_calib[i],input$parameter)
      if(input[row_num, 'distribution_type']=='gaussian') {
        par_calib_subsample[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      } else if(input[row_num, 'distribution_type']=='lognormal') {
        par_calib_subsample[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
      } else {
        print('ERROR - unknown prior distribution type')
      }
    }
    par_calib[ind_subsample,] <- par_calib_subsample

    # run model, precalibration windows
    model_out <- sapply(1:n_sample_max, function(ss) {
        model_forMCMC(par_calib=par_calib_subsample[ss,],
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
    for (ss in 1:n_sample_max) {
      if( any(model_out[,ss] < 100) |
          any(model_out[,ss] > 1e4) |
          model_out[58,ss] < 280 | model_out[58,ss] > 400) {
        ibad <- c(ibad, ss)
      }
    }
    good_parameters[[k+1]] <- par_calib_subsample[-ibad,]
  }

}

parameters_good <- NULL
for (k in 1:(n_full + n_part)) {
  parameters_good <- rbind(parameters_good, good_parameters[[k]])
}
colnames(parameters_good) <- parnames_calib

tend <- proc.time()
print(paste('LHS precalibration simulations and filtering took ',tend[3]-tbeg[3],' seconds', sep=''))

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
filename_out <- paste('../output/geocarb_precalibration_parameters_alpha',100*alpha,'_sens',sens,'_',today,'.csv', sep='')
write.csv(rbind(parameters_good, bandwidths), file=filename_out, row.names=FALSE)
##==============================================================================


##==============================================================================
## Function for sampling from KDEs (inflating a LHS?)
##===================================================

## Read KDE results file, separate into parameters and the bandwidths
filename_in <- filename_out
#alpha <- 0; filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL2_24Mar2018.csv'
#alpha <- 0.10; filename_in <- '../output/geocarb_precalibration_parameters_alpha10_sensL2_25Mar2018.csv'
#alpha <- 0.34; filename_in <- '../output/geocarb_precalibration_parameters_alpha34_sensL2_24Mar2018.csv'
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

## Get a reference simulation for integrated sensitivity measure
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


## Sobol' wrapper - assumes uniform distributions on parameters

source('GEOCARB_sensitivity_co2.R')

geocarb_sobol_co2_ser <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, ind_expected_time,
                          ind_expected_const, iteration_threshold) {
  co2_output <- sensitivity_co2(par_calib_scaled, l_scaled=TRUE,
                par_fixed=par_fixed, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input, ind_expected_time=ind_expected_time,
                ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold, model_ref=model_ref, sens=sens)
#  co2_output_centered <- co2_output - mean(co2_output, na.rm=TRUE)
  co2_output_centered <- co2_output - mean(co2_output[is.finite(co2_output)])
  co2_output_centered[is.infinite(co2_output)] <- 0
  return(co2_output_centered)
}

n_sample <- 500
n_bootstrap <- 5000
Ncore <- 1

## Sample parameters (need 2 data frames)
#parameters_sample1 <- kde_sample(n_sample, parameters_node, bandwidths)
#parameters_sample2 <- kde_sample(n_sample, parameters_node, bandwidths)
## Or take directly from the precalibration?
n_half <- floor(0.5*nrow(parameters_node))
parameters_sample1 <- parameters_node[1:n_half,]
parameters_sample2 <- parameters_node[(n_half+1):(2*n_half),]
colnames(parameters_sample1) <- colnames(parameters_sample2) <- parnames_calib

## Actually run the Sobol'
t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_co2_ser,
                           parameters_sample1,
                           parameters_sample2,
                           scheme='B',
                           nboot=n_bootstrap,
                           par_fixed=par_fixed0, parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                           ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                           input=input,
                           ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold))

print(paste('Sobol simulations took ',t.out[3],' seconds', sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.sobol <- paste('../output/sobol_alpha',100*alpha,'_',today,'.rds', sep='')
saveRDS(s.out, filename.sobol)

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('max. conf int=',max_conf_int, ' want less than:  0.1 * max. sensitivity index=',max_sens_ind, sep=''))
##==============================================================================



##==============================================================================
## Write indices, results we'd need to file

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_alpha',100*alpha,'_sens',sens,'_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_alpha',100*alpha,'_sens',sens,'_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames_calib,
                                     s.out$S[,1],
                                     s.out$S[,4],
                                     s.out$S[,5],
                                     s.out$T[,1],
                                     s.out$T[,4],
                                     s.out$T[,5]))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

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
##==============================================================================


##==============================================================================
## End
##==============================================================================
