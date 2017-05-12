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
  par,
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

  # Check lower/upper bounds first
  if( all(par <= bounds_calib[,'upper']) & all(par >= bounds_calib[,'lower']) ){

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
          lpri_new <- dnorm(x=par[ind_time_calib[((i-1)*ageN+1):(i*ageN)]], mean=time_arrays[,col_num], sd=(0.5*time_arrays[,col_num+1]), log=TRUE)
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          lpri_new <- dlnorm(x=par[ind_const_calib[((i-1)*ageN+1):(i*ageN)]], meanlog=log(time_arrays[,col_num]), sdlog=log(0.5*time_arrays[,col_num+1]), log=TRUE)
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
          lpri_new <- dnorm(x=par[ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          lpri_new <- dlnorm(x=par[ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else {
          print('ERROR - unknown prior distribution type')
        }
        lpri_const <- lpri_const + lpri_new
      }
    }

  } else {
    # if any of the calibration parameters are outside of the prior range
    lpri <- -Inf
  }

  # Add them up, assuming independence
  lpri <- lpri_time + lpri_const

  return(lpri)
}
##==============================================================================



##==============================================================================
## Likelihood function
##====================
log_like <- function(
  par,
  par_fixed,
  parnames_calib,
  parnames_fixed,
  age,
  ageN,
  ind_const_calib,
  ind_time_calib,
  ind_const_fixed,
  ind_time_fixed
){

  llike <- 0

  # run the model
  model_out <- model_forMCMC(par=par,
                             par_fixed=par_fixed,
                             parnames_calib=parnames_calib,
                             parnames_fixed=parnames_fixed,
                             age=age,
                             ageN=ageN,
                             ind_const_calib=ind_const_calib,
                             ind_time_calib=ind_time_calib,
                             ind_const_fixed=ind_const_fixed,
                             ind_time_fixed=ind_time_fixed)

  # compare against data

# TONY TODO
# TONY TODO
# TONY TODO

  return(llike)
}
##==============================================================================



##==============================================================================
## Posterior distribution
##=======================
log_post <- function(
  par,
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

  # calculate log-prior probability at these parameter values
  lpri <- log_prior(par=par,
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

# TONY TODO
# TONY TODO
# TONY TODO

  }

  # combine
  lpost <- lpri + llike

  return(lpost)
}
##==============================================================================



##==============================================================================
## End
##==============================================================================
