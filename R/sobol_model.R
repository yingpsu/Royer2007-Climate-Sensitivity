##==============================================================================
## sobol_model.R
##
## The Sobol integrals will all choke with NAs, which occur with
## some regularity due to the `failed_runs` in GEOCARB.  So, we just do the
## calculation ourselves for the first-order and total sensitivity indices.
##
## Assumes the parameters are already appropriately scale for use in the model.
##
##  sobol_model_par -- parallel version. first creates one large array of all of
##                     the parameter sets at which we need to run the model.
##                     this saves on the overhead of starting/stopping a cluster
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


##==============================================================================
sobol_model <- function(parameters, sens,
                        par_fixed, parnames_calib, parnames_fixed, parnames_time,
                        age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                        ind_time_fixed, ind_expected_time, ind_expected_const,
                        iteration_threshold, input, model_ref=NULL, data_calib=NULL,
                        do_sample_tvq=FALSE, par_time_center=NULL, par_time_stdev=NULL,
                        n_core=1, export_names=NULL
                        ){

  n_simulations <- nrow(parameters)

  if(is.null(n_simulations)) {
    # run the model ensemble under the sample `parameters`
    # Single simulation case
    model_out <- model_forMCMC(par_calib=parameters,
                               par_fixed=par_fixed,
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
                               do_sample_tvq=do_sample_tvq,
                               par_time_center=par_time_center,
                               par_time_stdev=par_time_stdev)[,'co2']

    # calculate sensitivity metric
    if (sens=='pres') {
      # present-day CO2
      model_sens <- model_out[ageN]
    } else if (sens=='L2') {
      # L2 norm
      model_sens <- sum((model_out-model_ref)^2)
    } else if (sens=='L1') {
      # L1 norm
      mod <- model_out[min(ind_mod2obs):max(ind_mod2obs)]
      icomp <- which(is.finite(mod) & !is.na(mod))
      ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
      model_sens <- sum(abs(mod[icomp]-ref[icomp]))
    } else if (sens=='NS') {
      # Nash-Sutcliffe efficiency
      model_stdy <- model_out[ind_mod2obs]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
      sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
      sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
      model_sens <- 1 - sse/sst
    } else if (sens=='NSL') {
      # Nash-Sutcliffe efficiency, on a logarithmic scale
      model_stdy <- model_out[ind_mod2obs]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy) & all(model_stdy > 0))  # only compare valid values
      sse <- sum( (log(model_stdy[icomp]) - log(data_calib$co2[icomp]))^2 )
      sst <- sum( (log(data_calib$co2[icomp]) - mean(log(data_calib$co2[icomp])))^2 )
      model_sens <- 1 - sse/sst
    }
  } else {
    # run the model ensemble under the sample `parameters`
    model_out <- sapply(1:n_simulations, function(ss) {
                   model_forMCMC(par_calib=parameters[ss,],
                                 par_fixed=par_fixed,
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
                                 do_sample_tvq=do_sample_tvq,
                                 par_time_center=par_time_center,
                                 par_time_stdev=par_time_stdev)[,'co2']})

    # calculate sensitivity metric
    if (sens=='pres') {
      # present-day CO2
      model_sens <- model_out[ageN,]
    } else if (sens=='L2') {
      # L2 norm
      model_sens <- apply(X=(model_out-model_ref)^2, MARGIN=2, FUN=sum)
    } else if (sens=='L1') {
      # L1 norm
      ##model_sens <- apply(X=abs(model_out[min(ind_mod2obs):max(ind_mod2obs),]-model_ref[min(ind_mod2obs):max(ind_mod2obs)]), MARGIN=2, FUN=sum)
      model_sens <- rep(NA, n_simulations)
      for (ss in 1:n_simulations) {
        mod <- model_out[min(ind_mod2obs):max(ind_mod2obs), ss]
        icomp <- which(is.finite(mod) & !is.na(mod))
        ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
        model_sens[ss] <- sum(abs(mod[icomp]-ref[icomp]))
      }
    } else if (sens=='NS') {
      # Nash-Sutcliffe efficiency
      model_sens <- rep(NA, n_simulations)
      for (ss in 1:n_simulations) {
        model_stdy <- model_out[ind_mod2obs, ss]
        icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
        sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
        sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
        model_sens[ss] <- 1 - sse/sst
      }
    } else if (sens=='NSL') {
      # Nash-Sutcliffe efficiency, on a logarithmic scale
      model_sens <- rep(NA, n_simulations)
      for (ss in 1:n_simulations) {
        model_stdy <- model_out[ind_mod2obs, ss]
        icomp <- which(is.finite(model_stdy) & !is.na(model_stdy) & all(model_stdy > 0))  # only compare valid values
        sse <- sum( (log(model_stdy[icomp]) - log(data_calib$co2[icomp]))^2 )
        sst <- sum( (log(data_calib$co2[icomp]) - mean(log(data_calib$co2[icomp])))^2 )
        model_sens[ss] <- 1 - sse/sst
      }
    }
  }

  return(model_sens)
}

##==============================================================================
sobol_model_par <- function(parameters, sens,
                            par_fixed, parnames_calib, parnames_fixed, parnames_time,
                            age, ageN, ind_const_calib, ind_time_calib,
                            ind_const_fixed, ind_time_fixed, ind_expected_time,
                            ind_expected_const, iteration_threshold, input,
                            model_ref=NULL, data_calib=NULL, do_sample_tvq=FALSE,
                            par_time_center=NULL, par_time_stdev=NULL,
                            n_core=1, export_names=NULL
                            ){

  n_simulations <- nrow(parameters)
  output <- rep(0, n_simulations)

  if(is.null(n_simulations)) {
    # run the model ensemble under the sample `parameters`
    # Single simulation case
    model_out <- model_forMCMC(par_calib=parameters,
                               par_fixed=par_fixed,
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
                               do_sample_tvq=do_sample_tvq,
                               par_time_center=par_time_center,
                               par_time_stdev=par_time_stdev)[,'co2']

    # calculate sensitivity metric
    if (sens=='pres') {
      # present-day CO2
      model_sens <- model_out[ageN]
    } else if (sens=='L2') {
      # L2 norm
      model_sens <- sum((model_out-model_ref)^2)
    } else if (sens=='L1') {
      # L1 norm
      mod <- model_out[min(ind_mod2obs):max(ind_mod2obs)]
      icomp <- which(is.finite(mod) & !is.na(mod))
      ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
      model_sens <- sum(abs(mod[icomp]-ref[icomp]))
    } else if (sens=='NS') {
      # Nash-Sutcliffe efficiency
      model_stdy <- model_out[ind_mod2obs]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
      sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
      sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
      model_sens <- 1 - sse/sst
    } else if (sens=='NSL') {
      # Nash-Sutcliffe efficiency, on a logarithmic scale
      model_stdy <- model_out[ind_mod2obs]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy) & all(model_stdy > 0))  # only compare valid values
      sse <- sum( (log(model_stdy[icomp]) - log(data_calib$co2[icomp]))^2 )
      sst <- sum( (log(data_calib$co2[icomp]) - mean(log(data_calib$co2[icomp])))^2 )
      model_sens <- 1 - sse/sst
    }

  } else {

    cores=detectCores()
    cl <- makeCluster(n_core)
    print(paste('Starting cluster with ',n_core,' cores', sep=''))
    registerDoParallel(cl)

    finalOutput <- foreach(i=1:n_simulations,
                           .combine=c,
                           .export=export_names,
                           .inorder=FALSE) %dopar% {

      dyn.load("../fortran/run_geocarb_unc.so")

      model_out <- model_forMCMC(par_calib=parameters[i,],
                                 par_fixed=par_fixed,
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
                                 do_sample_tvq=do_sample_tvq,
                                 par_time_center=par_time_center,
                                 par_time_stdev=par_time_stdev)[,'co2']

      # calculate sensitivity metric
      if (sens=='pres') {
        # present-day CO2
        model_sens <- model_out[ageN]
      } else if (sens=='L2') {
        # L2 norm
        model_sens <- sum((model_out-model_ref)^2)
      } else if (sens=='L1') {
        # L1 norm
        mod <- model_out[min(ind_mod2obs):max(ind_mod2obs)]
        icomp <- which(is.finite(mod) & !is.na(mod))
        ref <- model_ref[min(ind_mod2obs):max(ind_mod2obs)]
        model_sens <- sum(abs(mod[icomp]-ref[icomp]))
      } else if (sens=='NS') {
        # Nash-Sutcliffe efficiency
        model_stdy <- model_out[ind_mod2obs]
        icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
        sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
        sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
        model_sens <- 1 - sse/sst
      } else if (sens=='NSL') {
        # Nash-Sutcliffe efficiency, on a logarithmic scale
        model_stdy <- model_out[ind_mod2obs]
        icomp <- which(is.finite(model_stdy) & !is.na(model_stdy) & all(model_stdy > 0))  # only compare valid values
        sse <- sum( (log(model_stdy[icomp]) - log(data_calib$co2[icomp]))^2 )
        sst <- sum( (log(data_calib$co2[icomp]) - mean(log(data_calib$co2[icomp])))^2 )
        model_sens <- 1 - sse/sst
      }

      output[i] <- model_sens
    }
    print(paste(' ... done.'))
    stopCluster(cl)

    finalOutput[which(is.infinite(finalOutput))] <- NA

    return(finalOutput)
  }
}

##==============================================================================
## End
##==============================================================================
