##==============================================================================
## fit_likelihood_surface_mixture.R
##
## Fits a single Gaussian distribution to all of the data for each time slice.
##
## Yields:
##   loglikelihood_smoothed(modeled_co2, likelihood_fit, idx_data) (function)
##   likelihood_fit        (list of functions for KDE fits for each time step)
##   idx_data              (model CO2 indices where there are data to compare)
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

##==============================================================================
## Get a standard model simulation
##================================

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
##==============================================================================


##==============================================================================
## Get calibration data
## - will be ignoring the fitted distributions
##============================================

# later, we will sample from these and at each time step, fit a KDE to the
# distribution of all of the samples within each time step
upper_bound_co2 <- .upper_bound_co2
lower_bound_co2 <- .lower_bound_co2

if (data_choice=="PR2011") {
    filename.data <- "../input_data/CO2_Proxy_PR2011_calib_SN-co2_28Aug2019.csv"
    source('getData_PR2011.R')
} else if (data_choice=="F2017") {
  if (dist=='sn') {
    filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
  } else if (dist=='nm' | dist=='ln') {
    filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_NM-co2_25Sep2018.csv'
  } else {
    print(paste("ERROR: unknown distribution type",dist))
  }
  # Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
  data_to_assim <- cbind( c("paleosols" , TRUE),
                          c("alkenones" , TRUE),
                          c("stomata"   , TRUE),
                          c("boron"     , TRUE),
                          c("liverworts", TRUE) )
  source('getData.R')
  ind_data    <- which(data_to_assim[2,]==TRUE)
  n_data_sets <- length(ind_data)
  ind_assim   <- vector("list",n_data_sets)
  for (i in 1:n_data_sets) {
    ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
  }
  data_calib <- data_calib_all[unlist(ind_assim),]
}
##==============================================================================


##==============================================================================
## Fit mixture of skew-normals, each data point is a component
## weights for each component computed as in Park and Royer (2011):
##   --> (pCO2)avg,k = exp[sum(log(pCO2_i)/sig2_i)/sum(1/sig2_i)]
##   --> sig_i = log((pCO2)upper,i/(pCO2)i)
##   --> sig2avg,k = sum( (log(pCO2_i/(pCO2)avg,k))^2 )/(K-1)
##       note that sig_avg,k is the uncertainty in log(pCO2)
##   --> K = # data point in time slice k
##   --> if K=1, then (pCO2)avg,k = pCO2_i and sig2avg,k = sig2_i
## Weights here are the (1/sig^2)/sum(1/sig^2)
##================================================================

time <- model_out[,1]
n_time <- length(time)
dtime <- median(diff(time))
likelihood_fit <- vector('list', n_time)
x_co2 <- seq(from=.lower_bound_co2, to=.upper_bound_co2, by=1)

if (lhood_choice=="mixture") {
  for (tt in 1:n_time) {
      idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
      if (length(idx) > 0) {
          sigmas <- log(data_calib$co2_high[idx]/data_calib$co2[idx])
          means <- data_calib$co2[idx]
          wgts <- (1/(sigmas^2))/sum(1/(sigmas^2))
          wgts <- wgts/sum(wgts) # to be safe!
          # fit skew-normal mixture model
          f_co2 <- vector('list', length(idx))
          for (ii in 1:length(idx)) {
              #f_co2[[ii]] <- dlnorm(x_co2, meanlog=log(data_calib$co2[idx[ii]]), sdlog=sigmas[ii]) # log-normals from Park and Royer (2011)
              if (dist=='sn') {
                f_co2[[ii]] <- dsn(x_co2,xi=data_calib$xi_co2[idx[ii]], omega=data_calib$omega_co2[idx[ii]], alpha=data_calib$alpha_co2[idx[ii]])
              } else if (dist=='ln') {
                f_co2[[ii]] <- dlnorm(x_co2, logmean=log(data_calib$co2[idx[ii]]), sdlog=log(data_calib$sigma_co2[idx[ii]]))
              } else if (dist=='nm') {
                  f_co2[[ii]] <- dnorm(x_co2, mean=data_calib$mu_co2[idx[ii]], sd=data_calib$sigma_co2[idx[ii]])
              } else {
                print(paste("ERROR: unknown distribution type",dist))
              }
          }
          f_co2_mix <- wgts[1]*f_co2[[1]]
          if (length(idx) > 1) {
              for (ii in 2:length(idx)) {f_co2_mix <- f_co2_mix + wgts[ii]*f_co2[[ii]]}
          }
          # fit linear interpolation around KDE
          likelihood_fit[[tt]] <- approxfun(x_co2, f_co2_mix)
      }
  }
} else if (lhood_choice=="unimodal") {
  for (tt in 1:n_time) {
      idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
      if (length(idx) > 0) {
          sigmas <- log(data_calib$co2_high[idx]/data_calib$co2[idx])
          means <- data_calib$co2[idx]
          center <- exp( sum(log(means)/(sigmas^2)) / sum(1/(sigmas^2)) )
          if (length(idx) > 1) {
              stdev <- sqrt( sum( log(means/center)^2 )/(length(idx)-1) )
          } else if (length(idx)==1) {
              stdev <- log(data_calib$co2_high[idx]/data_calib$co2_low[idx])
          }
          # fit linear interpolation around the log-normal
          likelihood_fit[[tt]] <- approxfun(x_co2, dlnorm(x_co2, meanlog=log(center), sdlog=stdev))
      }
  }
} else {
  print(paste("ERROR: unknown lhood_choice",lhood_choice))
}
##==============================================================================


##==============================================================================
## Which model time steps have data?
##==================================
idx_data <- NULL
for (tt in 1:n_time) {
    if(!is.null(likelihood_fit[[tt]])) {idx_data <- c(idx_data,tt)}
}
n_data_fit <- length(idx_data)
##==============================================================================


##==============================================================================
## Function to evaluate the log-likelihood
##==============================================================================
loglikelihood_smoothed <- function(modeled_co2, likelihood_fit, idx_data, stdev=NULL) {
  llike <- 0
  if(!is.null(stdev)) {
    dx <- 10
    co2 <- seq(from=-10000,to=15000,by=dx) # excessive width to make sure nothing is clipped
    for (ii in idx_data) {
      f_like <- likelihood_fit[[ii]](co2); f_like[is.na(f_like)] <- 0
      f_unc <- dnorm(co2, mean=0, sd=stdev)
      f_conv <- convolve(f_like*dx, f_unc*dx)
      # re-order to center the convolution
      i0 <- which(co2==0)
      f_conv <- c(f_conv[(length(f_conv)-i0):length(f_conv)], f_conv[1:(length(f_conv)-i0-1)])
      if(is.null(ncol(modeled_co2))) { #print('1D')
        llike <- llike + log(approx(co2,f_conv,xout=modeled_co2[ii])$y)
      } else { #print("2D")
        llike <- llike + log(approx(co2,f_conv,xout=modeled_co2[ii,'co2'])$y)
      }
      if (is.na(llike) | is.infinite(llike)) {return(-Inf)}
    }
  } else {
    for (ii in idx_data) {
      llike <- llike + log(likelihood_fit[[ii]](modeled_co2[ii]))
      if (is.na(llike) | is.infinite(llike)) {return(-Inf)}
    }
  }
  return(llike)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
