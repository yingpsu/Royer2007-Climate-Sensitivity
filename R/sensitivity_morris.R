##==============================================================================
## sensitivity_morris.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
## Using Method of Morris one-at-a-time (sometimes referred to as MOAT)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

## Send in bounds of [0,1] for each parameter, and then scale within the
## `log_like_sensitivity` function to the parameters' distributions (using CDF)

alpha <- 0.1
perc_lower <- rep(0.5*alpha  , n_parameters)
perc_upper <- rep(1-0.5*alpha, n_parameters)

repetitions <- c(10, 20, 40, 80)
levels <- c(10, 20, 40, 80)

n_rep <- length(repetitions)
n_lev <- length(levels)

names_rep <- paste('r',repetitions[1],sep='')
if(n_rep>1) {for (rr in 2:n_rep) {names_rep <- c(names_rep, paste('r',repetitions[rr], sep=''))}}
names_lev <- paste('l',levels[1],sep='')
if(n_lev>1) {for (ll in 2:n_lev) {names_lev <- c(names_lev, paste('l',repetitions[ll], sep=''))}}

mom <- vector('list', n_rep); names(mom) <- names_rep
for (rr in names_rep) {
  mom[[rr]] <- vector('list', n_lev); names(mom[[rr]]) <- names_lev
}

## initialize other output
mu <- mu_star <- med <- med_star <- sigma <- stderr_mu <- mom_table <- names_signif <- mom

if(ldomorris) {
  for (rr in 1:n_rep) {
    for (ll in 1:n_lev) {
      tbeg <- proc.time()
      grid_jump <- round(0.5*levels[ll])
      mom[[rr]][[ll]] <- morris(model=log_like_sensitivity, factors=parnames_calib, r=repetitions[rr],
                                binf=perc_lower, bsup=perc_upper, scale=FALSE,
                                design=list(type="oat", levels=levels[ll], grid.jump=grid_jump),
                                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                                input=input, data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                                iteration_threshold=iteration_threshold)
      tend <- proc.time()
      print(paste(length(mom[[rr]][[ll]]$y),' simulations took ',(tend-tbeg)[3]/60,' minutes',sep=''))
      mu[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, mean)
      mu_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) mean(abs(x)))
      med[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, median)
      med_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) median(abs(x)))
      sigma[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, sd)
      stderr_mu[[rr]][[ll]] <- sigma[[rr]][[ll]]/sqrt(repetitions[rr])
      mom_table[[rr]][[ll]] <- cbind(mu[[rr]][[ll]], mu_star[[rr]][[ll]], sigma[[rr]][[ll]], stderr_mu[[rr]][[ll]])
    }
  }
  ## Get the set of significant (mu_star > 2*stderr_mu) parameters for each of the
  ## simulations; take the intersection (union?) of the "okay" ones as the
  ## parameter set for calibration?
  names_all <- NULL
  for (rr in 1:n_rep) {
    for (ll in 1:n_lev) {
      names_signif[[rr]][[ll]] <- parnames_calib[which(mu_star[[rr]][[ll]] > 2*stderr_mu[[rr]][[ll]])]
      names_all <- union(names_all, names_signif[[rr]][[ll]])
    }
  }
}
##==============================================================================
## End
##==============================================================================
