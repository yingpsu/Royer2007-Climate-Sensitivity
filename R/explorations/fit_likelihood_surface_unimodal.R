##==============================================================================
## fit_likelihood_surface_unimodal.R
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

##==============================================================================
## Get a standard model simulation
##================================

if(DO_SAMPLE_TVQ) {
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
} else {
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
                             iteration_threshold=iteration_threshold)
}
##==============================================================================


##==============================================================================
## Get calibration data
## - will be ignoring the fitted distributions
##============================================

# later, we will sample from these and at each time step, fit a KDE to the
# distribution of all of the samples within each time step
upper_bound_co2 <- .upper_bound_co2
lower_bound_co2 <- .lower_bound_co2

filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

source('GEOCARB-2014_getData.R')

ind_data    <- which(data_to_assim[2,]==TRUE)
n_data_sets <- length(ind_data)
ind_assim   <- vector("list",n_data_sets)
for (i in 1:n_data_sets) {
  ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
}

data_calib <- data_calib_all[unlist(ind_assim),]
##==============================================================================


##==============================================================================
## Fit normals for each time slice as done in Park and Royer (2011):
##   --> (pCO2)avg,k = exp[sum(log(pCO2_i)/sig2_i)/sum(1/sig2_i)]
##   --> sig_i = log((pCO2)upper,i/(pCO2)i)
##   --> sig2avg,k = sum( (log(pCO2_i/(pCO2)avg,k))^2 )/(K-1)
##       note that sig_avg,k is the uncertainty in log(pCO2)
##   --> K = # data point in time slice k
##   --> if K=1, then (pCO2)avg,k = pCO2_i and sig2avg,k = sig2_i
##================================================================

time <- model_out[,1]
n_time <- length(time)
dtime <- median(diff(time))
likelihood_fit <- vector('list', n_time)

for (tt in 1:n_time) {
    idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
    if (length(idx) > 0) {
        sigmas <- log(data_calib$co2_high[idx]/data_calib$co2[idx])
        means <- data_calib$co2[idx]
        center <- exp( sum(log(means)/(sigmas^2)) / sum(1/(sigmas^2)) )
        if (length(idx) > 1) {
            stdev <- sqrt( sum( log(means/center)^2 )/(length(idx)-1) )
        }
        # fit Gaussian with given mean=center, sd=stdev
        # hack-ish, but works well with how the mixture likelihood is computed
        x_co2 <- seq(from=lower_bound_co2, to=upper_bound_co2, by=1)
        ##Old and wrong:
        ##samples <- dnorm(log(co2), mean=log(center), sd=stdev)
        ##samples <- exp(samples) # since log(pCO2) assumed to be ~ normal
        # fit linear interpolation around KDE
        likelihood_fit[[tt]] <- approxfun(x_co2, dlnorm(x_co2, meanlog=log(center), sdlog=stdev))
    }
}

if(FALSE) {
# check out a cross-sectional slice
co2 <- seq(from=1,to=10000,by=10)

# a representative cross-sectional plot
t <- 34
f_fit <- likelihood_fit[[t]](co2)
plot(co2,f_fit, xlab='CO2 (ppmv)', ylab='density', type='l', xlim=c(0,6000))
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

if(FALSE){
##OLD
loglikelihood_smoothed <- function(modeled_co2, likelihood_fit, idx_data) {
  llike <- 0
  for (ii in idx_data) {
    llike <- llike + log(likelihood_fit[[ii]](modeled_co2[ii]))
    if (is.na(llike) | is.infinite(llike)) {return(-Inf)}
  }
  return(llike)
}
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
