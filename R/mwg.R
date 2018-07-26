##==============================================================================
## mwg.R
##
## Actual Metropolis-within-Gibbs sampling.
##
##  scale -- individual transition proposal variances
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

#mcmc_mwg <- function(log_post, n, init, scale, adapt, n.start, parallel, ncore,
#                     par_fixed, parnames_calib, parnames_fixed, age, ageN,
#                     ind_const_calib, ind_time_calib, ind_const_fixed,
#                     ind_time_fixed, input, time_arrays, bounds_calib,
#                     data_calib, ind_mod2obs, ind_expected_time,
#                     ind_expected_const, iteration_threshold) {

mcmc_mwg <- function(log_post, n, init,
                     scale, adapt=FALSE, n.start=0,
                     parallel=FALSE, n.core=1, ...) {

  # number of parameters to sample
  np <- length(init)

  # set up array for the samples and acceptances
  samples <- mat.or.vec(nr=n, nc=np)
  samples[1,] <- init
  n_accept <- rep(0, np)

  # get log-posterior score at the old parameters
  lp_old <- log_post(samples[1,], ...)

  # main iteration
  if(adapt) {
    # with adaptation of the proposal variances
    pb <- txtProgressBar(min=0,max=n,initial=0,style=3)
    for (t in 2:n) {
      # adaptation of variances; d=1 (Haario et al 2001) because each Metropolis
      # iteration is independent of the others
      if (t > n.start) {scale <- (2.4*2.4) * diag(cov(samples[1:(t-1),]))}
      # iterate over the parameters
      for (p in 1:np) {
        # new parameters are mostly the old ones...
        par_new <- samples[t-1,]
        # ... but replacing the pth one with a new step
        par_new[p] <- rnorm(mean=samples[t-1,p], sd=sqrt(scale[p]), n=1)
        # get log-posterior score at new parameters
        lp_new <- log_post(par_new, ...)
        # get log of the metropolis acceptance probability
        if (log(runif(1)) < lp_new - lp_old) {
          # accept!
          lp_old <- lp_new
          samples[t,p] <- par_new[p]
          n_accept[p] <- n_accept[p]+1
        } else {
          # reject!
          samples[t,p] <- samples[t-1,p]
        }
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  } else {
    # no adaptation of the proposal variances
    pb <- txtProgressBar(min=0,max=n,initial=0,style=3)
    for (t in 2:n) {
      # iterate over the parameters
      for (p in 1:np) {
        # new parameters are mostly the old ones...
        par_new <- samples[t-1,]
        # ... but replacing the pth one with a new step
        par_new[p] <- rnorm(mean=samples[t-1,p], sd=sqrt(scale[p]), n=1)
        # get log-posterior score at new parameters
        lp_new <- log_post(par_new, ...)
        # get log of the metropolis acceptance probability
        if (log(runif(1)) < lp_new - lp_old) {
          # accept!
          lp_old <- lp_new
          samples[t,p] <- par_new[p]
          n_accept[p] <- n_accept[p]+1
        } else {
          # reject!
          samples[t,p] <- samples[t-1,p]
        }
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  }

  # output
  mwg_out <- vector('list', 2)
  names(mwg_out) <- c('samples','p_accept')
  mwg_out[[1]] <- samples
  mwg_out[[2]] <- n_accept/n

  return(mwg_out)
}

##==============================================================================
## End
##==============================================================================
