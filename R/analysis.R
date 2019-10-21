##==============================================================================
## analysis.R
##
## These codes should be run after `process_results.R`, which will yield the 
## chains_analysis_[datestamp].rds file named below, and before `plotting.R`,
## which will generate the plots for the paper based on the analysis conducted
## in this routine.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')
library(sn)
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_analysis <- paste('../output/analysis_',today,'.RData', sep="")

# Get whatever we need in order to plot various things.
# Then, save all on one RData file to be read by plotting routine.
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

chains <- readRDS("../output/chains_analysis_05Oct2019.rds")
# normal distribution results
parameters_nm <- chains$dFpAUsRlMnm
# control results
parameters <- chains$dFpAUsRlMsn
n_ensemble <- nrow(parameters)
n_parameter <- ncol(parameters)

##==============================================================================

param_choice <- 'all_stdev'  # Calibrate all 69 parameters? ("all") or only the 6 from Park and Royer 2011 ("PR2011")
data_choice <- 'F2017'       # Which data set?  PR2011 = Park and Royer (2011), or F2017 = Foster et al (2017)
lhood_choice <- 'mixture'    # Mixture model ("mixture") or unimodal ("unimodal")?
fSR_choice <- 'DT2019'       # Which fSR time series? ("PR2011", "LENTON", "DT2019")
dist <- 'sn'                 # kernel choice for each data point (sn (skew-normal), ln (log-normal), nm (normal))

source("model_setup.R")

##==============================================================================



##==============================================================================
## Further analysis of chains, formerly in process_results.R

##
## compute ESS quantiles
##

quantiles_i_want <- c(.025, 0.05, 0.17, 0.25, 0.5, 0.75, 0.83, 0.95, 0.975)

ess <- vector('list', 2)
names(ess) <- c("ess", "ess_glac")
for (ee in names(ess)) {
  ess[[ee]] <- mat.or.vec(length(chains), length(quantiles_i_want))
  rownames(ess[[ee]]) <- names(chains)
  colnames(ess[[ee]]) <- quantiles_i_want
  for (experiment in names(chains)) {
    if (ncol(chains[[experiment]])==7) {
      idx_ess <- 5
    } else {
      idx_ess <- 10
    }
    if (ee=="ess") {
      ess[[ee]][experiment,] <- quantile(chains[[experiment]][,idx_ess], quantiles_i_want)
    } else {
      ess[[ee]][experiment,] <- quantile(chains[[experiment]][,idx_ess]*chains[[experiment]][,idx_ess+1], quantiles_i_want)
    }
  }
}

## Get the quantiles from Park and Royer 2011:
pr2011_dat <- read.csv('../input_data/ParkRoyer2011_Fig3_85varred.csv')
pr2011_cdf <- approxfun(pr2011_dat[,4], pr2011_dat[,1])
pr2011_pdf <- approxfun(pr2011_dat[,1], pr2011_dat[,3])
x_pr2011 <- pr2011_cdf(quantiles_i_want)
ess$ess <- rbind(x_pr2011, ess$ess)

##
## compute model quantiles
##

library(gdata)
dat <- read.xls("../input_data/GEOCARB_output--both error envelopes_TW.xlsx", sheet="GEOCARB_output_INNER_BANDS")
model_quantiles_r2014 <- cbind(dat[,"X0.025"], dat[,"CO2"], dat[,"X0.975"])
colnames(model_quantiles_r2014) <- c('0.025','0.5','0.975')

model_experiment_quantiles <- vector('list', length(chains))
names(model_experiment_quantiles) <- names(chains)
for (ee in names(model_experiment_quantiles)) {
  if(substr(ee, 2, 2)=="P") {data_choice <- "PR2011"
  } else {data_choice <- "F2017"}
  if(substr(ee, 4, 4)=="P") {param_choice <- "PR2011_stdev"
  } else {param_choice <- "all_stdev"}
  if(substr(ee, 7, 7)=="O") {fSR_choice <- "PR2011"
  } else if(substr(ee, 7, 7)=="L") {fSR_choice <- "LENTON"
  } else if(substr(ee, 7, 7)=="R") {fSR_choice <- "DT2019"}
  if(substr(ee, 9, 9)=="U") {lhood_choice <- "unimodal"
  } else {lhood_choice <- "mixture"}
  dist <- substr(ee, 10, 11)
  source('model_setup.R')
  if(data_choice=="PR2011") {data_calib_pr2011 <- data_calib
  } else if(data_choice=="F2017") {data_calib_f2017 <- data_calib}
  model_out <- sapply(X=1:nrow(chains[[ee]]),
                FUN=function(k){model_forMCMC(par_calib=chains[[ee]][k,],
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
  model_experiment_quantiles[[ee]] <- mat.or.vec(nr=n_time, nc=3)
  for (t in 1:n_time) {
    model_experiment_quantiles[[ee]][t,] <- quantile(model_out[t,], c(.025, .5, .975))
  }
  colnames(model_experiment_quantiles[[ee]]) <- c('0.025', '0.5', '0.975')
}
model_experiment_quantiles$r2014 <- model_quantiles_r2014

##==============================================================================



##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

# need likelihood/posterior functions, to get max. posterior score simulation
source('calib_likelihood_unc.R')

# run the ensemble
model_out <- sapply(X=1:n_ensemble,
              FUN=function(k){model_forMCMC(par_calib=parameters[k,],
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
n_time <- nrow(model_out)

# get 5-95% range and median  are cols 1-3; max-post will be 4
quantiles_i_want <- c(0,0.005,.025,.05,.5,.95,.975,0.995,1)
model_quantiles <- mat.or.vec(nr=n_time, nc=(length(quantiles_i_want)+1))
colnames(model_quantiles) <- c('q000','q005','q025','q05','q50','q95','q975','q995','q100','maxpost')
for (t in 1:n_time) {
    #model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_out_noisy[t,], quantiles_i_want)
    model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_out[t,], quantiles_i_want)
}

# get posterior scores

lpost_out <- sapply(X=1:n_ensemble,
              FUN=function(k){log_post(par_calib=parameters[k,],
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

# needed for plotting
idx_gastaldo <- which(data_calib_all$reference=="Gastaldo et al., 2014")

# get a reference model, uncalibrated

model_ref <- model_forMCMC(par_calib=par_calib0,
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
                           par_time_stdev=par_time_stdev)[,'co2']

##======================================

## read results from Royer et al 2014 (DOI: 10.2475/09.2014.01, AJS), to compare
## the width of the CO2 hindcast quantiles

library(gdata)

dat <- read.xls("../input_data/GEOCARB_output--both error envelopes_TW.xlsx", sheet="GEOCARB_output_INNER_BANDS")
model_quantiles_royer <- cbind(dat[,"X0.025"], dat[,"CO2"], dat[,"X0.975"])
colnames(model_quantiles_royer) <- c('q025','co2','q975')

##==============================================================================



##==============================================================================
# Figure 3. Posterior probability density for Earth system sensitivity parameter
# (deltaT2X), relative to previous studies.

# above, have parameters[,ics]
ics <- match('deltaT2X', parnames_calib)
deltaT2X_density <- density(parameters[,ics], from=0, to=10)

iglac <- match('GLAC', parnames_calib)
glac_density <- density(parameters[,iglac], from=1, to=5)

deltaT2Xglac_density <- density(parameters[,ics]*parameters[,iglac], from=0, to=20)


#plot(deltaT2X_density$x, deltaT2X_density$y, type='l')

# distributions of deltaT2X from other studies?

# Royer et al 2007:  1.5 and 6.2 deg C (5–95% range)
# Park and Royer 2011:
pr2011_dat <- read.csv('../input_data/ParkRoyer2011_Fig3_85varred.csv')
pr2011_cdf <- approxfun(pr2011_dat[,1], pr2011_dat[,4])
pr2011_icdf <- approxfun(pr2011_dat[,4], pr2011_dat[,1])
pr2011_pdf <- approxfun(pr2011_dat[,1], pr2011_dat[,3])

deltaT2X_density_pr2011 <- vector('list', 2); names(deltaT2X_density_pr2011) <- c('x','y')
deltaT2X_density_pr2011$x <- deltaT2X_density$x
deltaT2X_density_pr2011$y <- pr2011_pdf(deltaT2X_density_pr2011$x)

# Park and Royer 2011 have about 16% probability above deltaT2X = 6 deg C
print(1-pr2011_cdf(6))
print(1-pr2011_cdf(7))

##======================================



##==============================================================================
# Figure 4. Radial convergence diagrams for the sensitivity experiment. The
# “Sensitive” parameters are those identified as significant in the sub-ensemble
# correlation evaluation. Filled blue nodes represent first-order sensitivity
# indices; filled purple nodes represent total-order sensitivity indices; filled
# gray bars represent second-order sensitivity indices for the interaction
# between the parameter pair.

# in plotting_sobol.R

##======================================



##==============================================================================
## Needed numbers

print("Median, 5-95%, 2.5-97.5%, 16-84% CIs:")
print("==== deltaT2X ====")
print(quantile(parameters[,which(parnames_calib=="deltaT2X")], c(.5,.05,.95, .025,.975, .16,.84)))
print("==== GLAC ====")
print(quantile(parameters[,which(parnames_calib=="GLAC")], c(.5,.05,.95, .025,.975)))
print("==== GLAC*deltaT2X ====")
print(quantile(parameters[,which(parnames_calib=="GLAC")]*parameters[,which(parnames_calib=="deltaT2X")], c(.5,.05,.95, .025,.975, .16,.84)))
print("")

Pr_CS_above_6 <- length(which(parameters[,which(parnames_calib=="deltaT2X")] > 6)) / nrow(parameters)
print(paste("Pr(deltaT2X > 6 degC) = ",Pr_CS_above_6, sep=""))
print("")
Pr_CS_above_7 <- length(which(parameters[,which(parnames_calib=="deltaT2X")] > 7)) / nrow(parameters)
print(paste("Pr(deltaT2X > 7 degC) = ",Pr_CS_above_7, sep=""))
print("")

Pr_CS_below_2.5 <- length(which(parameters[,which(parnames_calib=="deltaT2X")] < 2.5)) / nrow(parameters)
print(paste("Pr(deltaT2X < 2.5 degC) = ",Pr_CS_below_2.5, sep=""))

print("")

print("And for the normal experiment...")
print("Median, 5-95%, 2.5-97.5% CIs:")
print("==== deltaT2X ====")
print(quantile(parameters_nm[,which(parnames_calib=="deltaT2X")], c(.5,.05,.95, .025,.975)))
print("==== GLAC ====")
print(quantile(parameters_nm[,which(parnames_calib=="GLAC")], c(.5,.05,.95, .025,.975)))
print("==== GLAC*deltaT2X ====")
print(quantile(parameters_nm[,which(parnames_calib=="GLAC")]*parameters_nm[,which(parnames_calib=="deltaT2X")], c(.5,.05,.95, .025,.975)))
print("")

##==============================================================================



# Supplementary Figures



##==============================================================================
# Figure S2.  Likelihood slice from multimodal period.

mm_example240 <- vector('list', 3)
names(mm_example240) <- c('co2','fit','age')
idx <- which(age==240)
mm_example240$age <- age[idx]
mm_example240$co2 <- seq(from=1,to=10000,by=1)
mm_example240$fit <- likelihood_fit[[idx]](mm_example240$co2)


mm_example50 <- vector('list', 3)
names(mm_example50) <- c('co2','fit','age')
idx <- which(age==50)
mm_example50$age <- age[idx]
mm_example50$co2 <- seq(from=1,to=10000,by=1)
mm_example50$fit <- likelihood_fit[[idx]](mm_example50$co2)


#plot(mm_example240$co2, mm_example240$fit, type='l', xlab='CO2 (ppmv)', ylab='density')


## More investigation of multi-modality

# check model ensemble during the 240 Mya time slice
fit <- density(model_out[34,])
mm_modeled <- mm_example240
mm_modeled$fit <- fit$y
mm_modeled$co2 <- fit$x


# check precalibration ensemble during this time slice
parameters_precal <- readRDS('../output/precal_parameters_25Jun2019.rds')
model_out_precal <- sapply(X=1:nrow(parameters_precal),
              FUN=function(k){model_forMCMC(par_calib=parameters_precal[k,],
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
fit2 <- density(model_out_precal[34,])
mm_precal <- mm_example240
mm_precal$co2 <- fit2$x
mm_precal$fit <- fit2$y


# check parameters straight from the priors
library(lhs)
parameters_lhs0 <- randomLHS(100000, 69)

## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
parameters_lhs <- parameters_lhs0  # initialize
colnames(parameters_lhs) <- parnames_calib
for (i in 1:n_const_calib) {
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    parameters_lhs[,i] <- qnorm(p=parameters_lhs0[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    parameters_lhs[,i] <- qlnorm(p=parameters_lhs0[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}
for (i in (n_const_calib+1):length(parnames_calib)) {
  parameters_lhs[,i] <- qbeta(p=parameters_lhs0[,i], shape1=5, shape2=5)
}

model_out_priors <- sapply(1:nrow(parameters_lhs), function(ss) {
                    model_forMCMC(par_calib=parameters_lhs[ss,],
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
fit3 <- density(model_out_priors[34,which(model_out_priors[34,] < 1e4)])
mm_priors <- mm_example240
mm_priors$co2 <- fit3$x
mm_priors$fit <- fit3$y


# check parameters only GYM and deltaT2X varying - explain the multimodality?
parameters_essgym <- parameters_lhs
for (pp in 1:length(parnames_calib)) {
  if (parnames_calib[pp]!='deltaT2X' & parnames_calib[pp]!='GYM') {
    parameters_essgym[,pp] <- par_calib0[pp]
  }
}

model_out_essgym <- sapply(1:nrow(parameters_essgym), function(ss) {
                    model_forMCMC(par_calib=parameters_essgym[ss,],
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
fit4 <- density(model_out_essgym[34,which(model_out_essgym[34,] < 1e4)])
mm_essgym <- mm_example240
mm_essgym$co2 <- fit4$x
mm_essgym$fit <- fit4$y


# comparison
# figure not used in paper but potentially useful to pick apart constraints
pdf(paste('../figures/multimodality_SOM.pdf',sep=''),width=4,height=3, colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.8,.2,.2))
plot(mm_modeled$co2, mm_modeled$fit, type='l', lwd=2, lty=3, xlim=c(0,7000), ylim=c(0,0.0009),
     xlab='CO2 (ppmv)', ylab='Probability density')
lines(mm_example240$co2, mm_example240$fit, lwd=2)
lines(mm_precal$co2, mm_precal$fit, lwd=2, lty=4)
lines(mm_priors$co2, mm_priors$fit, lwd=2, lty=2)
#lines(mm_essgym$co2, mm_essgym$fit, lwd=2, lty=2, col='purple')
legend(3500, 0.0009, c('Likelihood','Priors','Precal','Posterior'),
       lty=c(1,2,4,3), col=c('black','black','black','black'), lwd=2, bty='n')
dev.off()

##======================================



##==============================================================================
# Figure for SOM?  Posterior model ensemble (gray shaded region denotes 5-95%
# credible range), maximum posterior score simulation (solid bold line) and
# uncalibrated model simulation (dashed line), with proxy data points
# superimposed (+ markers), assuming a symmetric (Gaussian) error structure for
# the proxy data as opposed to skew-normal (main text).

# run the ensemble
model_out_nm <- sapply(X=1:nrow(parameters_nm),
              FUN=function(k){model_forMCMC(par_calib=parameters_nm[k,],
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

# get 5-95% range and median  are cols 1-3; max-post will be 4
model_quantiles_nm <- mat.or.vec(nr=n_time, nc=(length(quantiles_i_want)+1))
colnames(model_quantiles_nm) <- colnames(model_quantiles)
for (t in 1:n_time) {
  model_quantiles_nm[t,1:length(quantiles_i_want)] <- quantile(model_out_nm[t,], quantiles_i_want)
}

# get posterior scores

lpost_out_nm <- sapply(X=1:nrow(parameters_nm),
              FUN=function(k){log_post(par_calib=parameters_nm[k,],
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

model_quantiles_nm[,'maxpost'] <- model_out_nm[,which.max(lpost_out_nm)]

##======================================



##==============================================================================
# Figure for SOM?  Posterior probability density for Earth system sensitivity
# parameter (deltaT2X), relative to previous studies), assuming a symmetric
# (Gaussian) error structure for the proxy data as opposed to skew-normal.

deltaT2X_density_nm <- density(parameters_nm[,ics], from=0, to=10)

#plot(deltaT2X_density$x, deltaT2X_density$y, type='l', lwd=2); lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y, type='l', lty=2, lwd=2)
##==============================================================================



##==============================================================================
# SOM figure -- likelihood surface and model quantiles
#===========

source('likelihood_surface_quantiles.R')

##==============================================================================



##==============================================================================
# SOM figure -- quantiles of deltaT2X as sample size increases
#===========

sample_sizes <- seq(from=250, to=nrow(chains$dFpAUsRlMsn), by=250)
sample_quantiles <- mat.or.vec(length(sample_sizes), 3)
for (i in 1:length(sample_sizes)) {
  sample_quantiles[i,] <- quantile(chains$dFpAUsRlMsn[1:sample_sizes[i], 10], c(.05,.5,.95))
}

##==============================================================================



##======================================
save.image(file=filename_analysis)
##======================================



##==============================================================================
## End
##==============================================================================
