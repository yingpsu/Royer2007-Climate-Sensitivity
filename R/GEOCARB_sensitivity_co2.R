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
  data_calib=NULL,
  sens
){

#saveRDS(par_calib_scaled, 'NS_par_debug.rds')

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

  if (is.null(n_simulations)) {n_simulations <- 1}

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
      model_sens <- model_out[ageN,]
    } else if (sens=='L2') {
      # L2 norm
      model_sens <- apply(X=(model_out-model_ref)^2, MARGIN=2, FUN=sum)
    } else if (sens=='L1') {
      # L1 norm
      ##model_sens <- apply(X=abs(model_out[min(ind_mod2obs):max(ind_mod2obs),]-model_ref[min(ind_mod2obs):max(ind_mod2obs)]), MARGIN=2, FUN=sum)
      model_sens <- rep(-999, n_simulations)
      for (ss in 1:n_simulations) {
        mod <- model_out[min(ind_mod2obs):max(ind_mod2obs), ss]
        icomp <- which(is.finite(mod) & !is.na(mod))
        ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
        model_sens[ss] <- sum(abs(mod[icomp]-ref[icomp]))
      }
    } else if (sens=='NS') {
      # Nash-Sutcliffe efficiency
      model_sens <- rep(-999, n_simulations)
      for (ss in 1:n_simulations) {
        model_stdy <- model_out[ind_mod2obs, ss]
        icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
        sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
        sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
        model_sens[ss] <- 1 - sse/sst

        # DEBUG
        #print(paste(ss,model_sens[ss],sse,sst,length(icomp)))
        if (is.infinite(model_sens[ss]) | is.na(model_sens[ss]) | sst==0 | length(icomp)==0) {print(paste('ss',ss,'-- sst',sst,'-- sse',sse,'-- model_sens',model_sens[ss],'-- model_out',model_out[,ss]))}
      }
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
      model_sens <- model_out[ageN]
    } else if (sens=='L2') {
      # L2 norm
      model_sens <- sum((model_out-model_ref)^2)
    } else if (sens=='L1') {
      # L1 norm
      ##model_sens <- sum(abs(model_out-model_ref))
      mod <- model_out[min(ind_mod2obs):max(ind_mod2obs)]
      icomp <- which(is.finite(mod) & !is.na(mod))  # only compare valid values
      ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
      model_sens <- sum(abs(mod[icomp]-ref[icomp]))
    } else if (sens=='NS') {
      # Nash-Sutcliffe efficiency
      model_stdy <- model_out[ind_mod2obs]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
      sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
      sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
      model_sens <- 1 - sse/sst

      # DEBUG
      #print(paste(ss,model_sens[ss],sse,sst,length(icomp)))
      if (is.infinite(model_sens) | is.na(model_sens) | sst==0) {print(paste('sst',sst,'-- sse',sse,'-- model_sens',model_sens,'-- model_out',model_out))}
    }
  }

  #ind_na <- which(is.na(model_sens))
  #if (length(ind_na)>0) {model_sens[ind_na] <- -Inf}
  #ind_inf <- which(is.infinite(model_sens))
  #if (length(ind_inf)>0) {model_sens[ind_inf] <- NA}
  #ind_low <- which(model_sens < -9)
  #if (length(ind_low)>0) {model_sens[ind_low] <- NA}

#  write.table(model_sens    , file='NS_sens_debug.txt', append=FALSE , sep = " ",
#              quote=FALSE    , row.names = FALSE , col.names=FALSE)

#print('debug here now')

  return(model_sens)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
