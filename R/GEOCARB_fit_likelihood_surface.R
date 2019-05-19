##==============================================================================
## GEOCARB_fit_likelihood_surface.R
##
## Read previously fit skew-normal distributions for each data point, and sample
## from all of them. Then, fit KDEs to each time step's data points, and build a
## linear interpolation for each of these (faster execution within MCMC).
##
## Requires:
##   dist                  one of: 'ga' (gamma), 'be' (beta), 'ln' (log-normal)
##                         or 'sn' (skew-normal), describing the distributions
##                         that were fit to each data point and its uncertain
##                         range in the processing step.
## Yields:
##   loglikelihood_smoothed(modeled_co2, likelihood_fit, idx_data) (function)
##   likelihood_fit        (list of functions for KDE fits for each time step)
##   idx_data              (model CO2 indices where there are data to compare)
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
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
## Get calibration data and fitted
## distributions to *each data point*
##===================================

# later, we will sample from these and at each time step, fit a KDE to the
# distribution of all of the samples within each time step
upper_bound_co2 <- .upper_bound_co2
lower_bound_co2 <- .lower_bound_co2

if(dist=='ga') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_GAMMA-co2_31Jul2018.csv'}
if(dist=='be') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_BETA-co2_13Sep2018.csv'}
if(dist=='ln') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_LN-co2_31Jul2018.csv'}
if(dist=='sn') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'}
if(dist=='nm') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_NM-co2_25Sep2018.csv'}
if(dist=='sn-100min') {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_100min_22Oct2018.csv'}
if(dist=='sn-mmrem')  {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_mmrem_27Oct2018.csv'}
if(dist=='nm-unifUnc')  {filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_NM-co2-unifUnc_29Nov2018.csv'}

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
## Sample and fit KDEs to each time step
##======================================

time <- model_out[,1]
n_time <- length(time)
dtime <- median(diff(time))
n_sample_per_point <- 10000
likelihood_fit <- vector('list', n_time)

for (tt in 1:n_time) {
    idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
    if (length(idx) > 0) {
        # sample from the distributions fit to each of the data points
        samples <- NULL
        for (ii in idx) {
            if (dist=='ga') {
                new_samples <- rgamma(shape=data_calib$shape_co2[ii], scale=data_calib$scale_co2[ii],
                                      n=n_sample_per_point)
            } else if (dist=='be') {
                new_samples <- rbeta(shape1=data_calib$shape1_co2[ii], shape2=data_calib$shape2_co2[ii],
                                     n=n_sample_per_point)
                new_samples <- lower_bound_co2+(upper_bound_co2-lower_bound_co2)*new_samples
            } else if (dist=='ln') {
                new_samples <- rlnorm(meanlog=data_calib$meanlog_co2[ii], sdlog=data_calib$sdlog_co2[ii],
                                      n=n_sample_per_point)
            } else if (dist=='sn' | dist=='sn-100min' | dist=='sn-mmrem') {
                new_samples <- rsn(xi=data_calib$xi_co2[ii], omega=data_calib$omega_co2[ii],
                                   alpha=data_calib$alpha_co2[ii], n=n_sample_per_point)
            } else if (dist=='nm' | dist=='nm-unifUnc') {
                new_samples <- rnorm(mean=data_calib$mu_co2[ii], sd=data_calib$sigma_co2[ii], n=n_sample_per_point)
            }
            samples <- c(samples, new_samples)
        }
        idx_filter <- which(samples < lower_bound_co2)
        if(length(idx_filter) > 0) {samples <- samples[-idx_filter]}
        # fit KDE
        density_fit <- density(samples, from=lower_bound_co2, to=upper_bound_co2)
        # fit linear interpolation around KDE
        likelihood_fit[[tt]] <- approxfun(density_fit)
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
  dx <- 10
  co2 <- seq(from=-10000,to=15000,by=dx) # excessive width to make sure nothing is clipped
  llike <- 0
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
