##==============================================================================
## GEOCARB_sensitivity_co2.R
##
## Todo...
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================


##==============================================================================
## function for sensitivity analysis
## (scales from [0,1] to the parameters' distributions)
##=====================================================
sensitivity_co2 <- function(
  par_calib_scaled, l_scaled=TRUE,
  par_fixed,
  parnames_calib,
  parnames_fixed,
  age,
  ageN,
  ind_const_calib,
  ind_time_calib,
  ind_const_fixed,
  ind_time_fixed,
  ind_expected_time,
  ind_expected_const,
  iteration_threshold,
  input,
  model_ref=NULL,
  sens
){

  # initialize
  par_calib <- par_calib_scaled

  # scale the parameters from [0,1] to their prior distributions IF needed
  n_const_calib <- length(ind_const_calib)
  n_simulations <- nrow(par_calib_scaled)

if(!l_scaled) {
  if (n_const_calib > 0) {
    if(n_simulations > 1) {
      for (i in 1:n_const_calib) {
        row_num <- match(parnames_calib[i],input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          par_calib[,i] <- qnorm(p=par_calib_scaled[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          par_calib[,i] <- qlnorm(p=par_calib_scaled[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
        } else {
          print('ERROR - unknown prior distribution type')
        }
      }
    } else {
      for (i in 1:n_const_calib) {
        row_num <- match(parnames_calib[i],input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          par_calib[i] <- qnorm(p=par_calib_scaled[ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          par_calib[i] <- qlnorm(p=par_calib_scaled[ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
        } else {
          print('ERROR - unknown prior distribution type')
        }
      }
    }
  }
}

  llike <- rep(0, n_simulations)

  # run the model
  if (n_simulations > 1) {
    model_out <- sapply(1:n_simulations, function(ss) {
                   model_forMCMC(par_calib=par_calib[ss,],
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
                                 iteration_threshold=iteration_threshold)[,'co2']})
    if (sens=='pres') {
      # present-day CO2
      model_present <- model_out[ageN,]
    } else if (sens=='L2') {
      # L2 norm
      model_present <- apply(X=(model_out-model_ref)^2, MARGIN=2, FUN=sum)
    } else if (sens=='L1') {
      # L1 norm
      model_present <- apply(X=abs(model_out-model_ref), MARGIN=2, FUN=sum)
    }
  } else {
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
    if (sens=='pres') {
      # present-day CO2
      model_present <- model_out[ageN]
    } else if (sens=='L2') {
      # L2 norm
      model_present <- sum((model_out-model_ref)^2)
    } else if (sens=='L1') {
      # L1 norm
      ##model_present <- sum(abs(model_out-model_ref))
      model_present <- sum(abs(model_out[min(ind_mod2obs):max(ind_mod2obs)]-model_ref[min(ind_mod2obs):max(ind_mod2obs)]))
    }
  }

  ind_na <- which(is.na(model_present))
  #if (length(ind_na)>0) {model_present[ind_na] <- -Inf}

  return(model_present)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
