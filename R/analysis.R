##==============================================================================
## analysis.R
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

rm(list=ls())

library(sn)
library(ncdf4)

# Get whatever we need in order to plot various things.
# Then, save all on one RData file to be read by plotting routine.



##==============================================================================
# Figure 1. Observations and fitted likelihood surface.

dist <- 'sn'
DO_SAMPLE_TVQ <- TRUE

data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_tvq_all.csv'

source('GEOCARB-2014_parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')
source('run_geocarbF.R')
source('GEOCARB_fit_likelihood_surface.R')

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

save.image(file='../output/analysis.RData')
##======================================



##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

ncdata <- nc_open('../output/geocarb_calibratedParameters_tvq_all_25Sep2018.nc')
parameters <- t(ncvar_get(ncdata, 'geocarb_parameters'))
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
n_ensemble <- nrow(parameters)
n_parameter <- ncol(parameters)

# need likelihood/posterior functions, to get max. posterior score simulation
source('GEOCARB-2014_calib_likelihood.R')

# run the ensemble
model_out <- sapply(X=1:n_ensemble,
              FUN=function(k){model_forMCMC(par_calib=parameters[k,],
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
                                            iteration_threshold=iteration_threshold,
                                            do_sample_tvq=DO_SAMPLE_TVQ,
                                            par_time_center=par_time_center,
                                            par_time_stdev=par_time_stdev)[,'co2']})
n_time <- nrow(model_out)

# get 5-95% range and median  are cols 1-3; max-post will be 4
model_quantiles <- mat.or.vec(nr=n_time, nc=4)
quantiles_i_want <- c(.05,.5,.95)
colnames(model_quantiles) <- c('q05','q50','q95','maxpost')
for (t in 1:n_time) {
  model_quantiles[t,1:3] <- quantile(model_out[t,], quantiles_i_want)
}

# get posterior scores

lpost_out <- sapply(X=1:n_ensemble,
              FUN=function(k){log_post(par_calib=parameters[k,],
                                      par_fixed=par_fixed0,
                                      parnames_calib=parnames_calib,
                                      parnames_fixed=parnames_fixed,
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

model_quantiles[,'maxpost'] <- model_out[,which.max(lpost_out)]

save.image(file='../output/analysis.RData')
##======================================



##==============================================================================
# Figure 3. Posterior probability density for Earth system sensitivity parameter
# (deltaT2X), relative to previous studies.

# above, have parameters[,ics]
ics <- match('deltaT2X', parnames)
deltaT2X_density <- density(parameters[,ics], from=0, to=10)

#plot(deltaT2X_density$x, deltaT2X_density$y, type='l')

# distributions of deltaT2X from other studies?

# TODO
# TODO
# TODO
# TODO

save.image(file='../output/analysis.RData')
##======================================



##==============================================================================
# Figure 4. Radial convergence diagrams for the sensitivity experiment. The
# “Sensitive” parameters are those identified as significant in the sub-ensemble
# correlation evaluation. Filled blue nodes represent first-order sensitivity
# indices; filled purple nodes represent total-order sensitivity indices; filled
# gray bars represent second-order sensitivity indices for the interaction
# between the parameter pair.

##======================================



# Supplementary Figures

##==============================================================================
# Figure S1.  Histograms and parameter pairs plots for the calibrated parameters.
# Use contours of probability density for the pairs plots.

##======================================



##==============================================================================
# Figure S2.  Posterior model ensemble (gray shaded region denotes 5-95%
# credible range), maximum posterior score simulation (solid bold line) and
# uncalibrated model simulation (dashed line), with proxy data points
# superimposed (+ markers), assuming a symmetric (Gaussian) error structure for
# the proxy data as opposed to skew-normal (main text).

##======================================



##==============================================================================
# Figure S3. Posterior probability density for Earth system sensitivity
# parameter (deltaT2X), relative to previous studies), assuming a symmetric
# (Gaussian) error structure for the proxy data as opposed to skew-normal.

##======================================



##==============================================================================
## End
##==============================================================================
