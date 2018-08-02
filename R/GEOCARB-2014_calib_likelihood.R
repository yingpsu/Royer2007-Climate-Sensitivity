##==============================================================================
## GEOCARB-2014_calib_likelihood.R
##
## Priors, likelihood function, posterior distribution
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================



##==============================================================================
## Prior distributions
##====================
log_prior <- function(
  par_calib,
  par_fixed,
  parnames_calib,
  parnames_fixed,
  age,
  ageN,
  ind_const_calib,
  ind_time_calib,
  ind_const_fixed,
  ind_time_fixed,
  input,
  time_arrays,
  bounds_calib
){

  lpri <- 0

  # Check lower/upper bounds first
  if( all(par_calib <= bounds_calib[,'upper']) & all(par_calib >= bounds_calib[,'lower']) ){

    # Gaussian process priors for the time-varying parameters
    # Might make simplifying assumption of independence between time slices
    lpri_time <- 0
    n_time_calib <- length(ind_time_calib)/ageN
    if (n_time_calib>0) {
      for (i in 1:n_time_calib) {
        lpri_new <- 0
        name <- parnames_calib[ind_time_calib[ageN*i]]
        row_num <- match(name,input$parameter)
        col_num <- match(name,colnames(time_arrays))
        if(input[row_num, 'distribution_type']=='gaussian') {
          lpri_new <- dnorm(x=par_calib[ind_time_calib[((i-1)*ageN+1):(i*ageN)]], mean=time_arrays[,col_num], sd=(0.5*time_arrays[,col_num+1]), log=TRUE)
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          lpri_new <- dlnorm(x=par_calib[ind_const_calib[((i-1)*ageN+1):(i*ageN)]], meanlog=log(time_arrays[,col_num]), sdlog=log(0.5*time_arrays[,col_num+1]), log=TRUE)
        } else {
          print('ERROR - unknown prior distribution type')
        }
        lpri_time <- lpri_time + sum(lpri_new)
      }
    }

    # Priors for the time-constant parameters
    lpri_const <- 0
    n_const_calib <- length(ind_const_calib)
    if (n_const_calib>0) {
      for (i in 1:n_const_calib) {
        row_num <- match(parnames_calib[i],input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          lpri_new <- dnorm(x=par_calib[ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          lpri_new <- dlnorm(x=par_calib[ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else {
          print('ERROR - unknown prior distribution type')
        }
        lpri_const <- lpri_const + lpri_new
      }
    }
    # Add them up, assuming independence
    lpri <- lpri_time + lpri_const
  } else {
    # if any of the calibration parameters are outside of the prior range
    lpri <- -Inf
  }

  return(lpri)
}
##==============================================================================



##==============================================================================
## Likelihood function
##====================
log_like <- function(
  par_calib,
  par_fixed,
  parnames_calib,
  parnames_fixed,
  age,
  ageN,
  ind_const_calib,
  ind_time_calib,
  ind_const_fixed,
  ind_time_fixed,
  data_calib,
  ind_mod2obs,
  ind_expected_time,
  ind_expected_const,
  iteration_threshold
){

  llike <- 0

  # run the model
  model_out <- model_forMCMC(par_calib=par_calib,
                             par_fixed=par_fixed,
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

  # use the same checks as the precalibration to get rid of unphysical simulations
  if(any(is.infinite(model_out)) | any(model_out < 0) | any(model_out > 1e4)) {
    llike <- -Inf
  } else {
    # compare against data
    # assumption of steady state in-between model time steps
    # note that these are not necessarily sequential in time
    model_stdy <- model_out[ind_mod2obs]
    llike <- sum( sapply(1:length(model_stdy), function(i) dsn(x=model_stdy[i],
                         xi=data_calib$xi_co2[i], omega=data_calib$omega_co2[i],
                         alpha=data_calib$alpha_co2[i], log=TRUE)) )
#    llike <- sum( sapply(1:length(model_stdy), function(i) dgamma(x=model_stdy[i],
#                         shape=data_calib$shape_co2[i], scale=data_calib$scale_co2[i],
#                         log=TRUE)) )
#    llike <- sum( sapply(1:length(model_stdy), function(i) dlnorm(x=model_stdy[i],
#                         meanlog=data_calib$meanlog_co2[i], sdlog=data_calib$sdlog_co2[i],
#                         log=TRUE)) )
    if(is.na(llike)) {llike <- -Inf}
  }
  return(llike)
}
##==============================================================================


##==============================================================================
## Negative posterior
##===================
neg_log_post <- function(par_calib, par_fixed, parnames_calib, parnames_fixed,
                         age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                         ind_time_fixed, input, time_arrays, bounds_calib, data_calib,
                         ind_mod2obs, ind_expected_time, ind_expected_const, iteration_threshold){

  lpost <- log_post(par_calib=par_calib, par_fixed=par_fixed, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
                data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold)
  return(-lpost)
}
##==============================================================================


##==============================================================================
## Posterior distribution
##=======================
log_post <- function(
  par_calib,
  par_fixed,
  parnames_calib,
  parnames_fixed,
  age,
  ageN,
  ind_const_calib,
  ind_time_calib,
  ind_const_fixed,
  ind_time_fixed,
  input,
  time_arrays,
  bounds_calib,
  data_calib,
  ind_mod2obs,
  ind_expected_time,
  ind_expected_const,
  iteration_threshold,
  n_shard=1
){

  lpri <- 0
  llike <- 0

  # calculate log-prior probability at these parameter values
  lpri <- log_prior(par_calib=par_calib,
                    par_fixed=par_fixed,
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
                    bounds_calib=bounds_calib)

  # calculate log-likelihood at these parameter values (with these data)
  # save time by not running model if prior is -Inf (i.e., parameters outside
  # prior ranges)
  # Fun note: it is faster to check is.finite() than ~is.infinite
  if(is.finite(lpri)) {
    llike <- log_like(par_calib=par_calib,
                      par_fixed=par_fixed,
                      parnames_calib=parnames_calib,
                      parnames_fixed=parnames_fixed,
                      age=age,
                      ageN=ageN,
                      ind_const_calib=ind_const_calib,
                      ind_time_calib=ind_time_calib,
                      ind_const_fixed=ind_const_fixed,
                      ind_time_fixed=ind_time_fixed,
                      data_calib=data_calib,
                      ind_mod2obs=ind_mod2obs,
                      ind_expected_time=ind_expected_time,
                      ind_expected_const=ind_expected_const,
                      iteration_threshold=iteration_threshold)
  }

  # combine
  lpost <- (1/n_shard)*lpri + llike
  #lpost <- lpri + llike

  return(lpost)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
