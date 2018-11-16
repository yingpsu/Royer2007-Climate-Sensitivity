#===============================================================================
#
#===============================================================================

setwd('~/codes/GEOCARB/output')
load('GEOCARB_MCMC_tvq_split_16Nov2018sn-split.RData')

# histograms figure
par(mfrow=c(5,4), mai=c(.3,.35,.1,.01))
for (ss in 1:3) {for (mm in 1:4) {hist(amcmc_par[[ss]][[mm]]$samples[5e5:1e6,10], freq=FALSE, xlim=c(1,9), main='', breaks=seq(0,10,by=0.25))}}

# history plots figure
par(mfrow=c(5,4), mai=c(.3,.35,.1,.01))
for (ss in 1:3) {for (mm in 1:4) {plot(amcmc_par[[ss]][[mm]]$samples[,10], type='l')}}

# grab posterior parameter ensembles from each shard
parameters <- vector('list', n_shard)
for (ss in 1:n_shard) {
    parameters[[ss]] <- amcmc_par[[ss]][[1]]$samples[seq((niter_mcmc*0.6+1),niter_mcmc,by=100),]
    for (mm in 2:length(amcmc_par[[1]])) {
        parameters[[ss]] <- rbind(parameters[[ss]],
                                  amcmc_par[[ss]][[mm]]$samples[seq((niter_mcmc*0.6+1),niter_mcmc,by=100),])
    }
}

# run posterior model ensembles for each shard
source('../R/run_geocarbF.R')
model_out <- vector('list',n_shard)
for (ss in 1:n_shard) {
    model_out[[ss]] <- sapply(1:nrow(parameters[[ss]]), function(ii) {
        model_forMCMC(par_calib=parameters[[ss]][ii,],
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
                      par_time_stdev=par_time_stdev)[,'co2']})
}

# posterior time slice at t=34 (age=240 Myr) from each shard
fits <- vector('list',n_shard)
tslice <- 34
for (ss in 1:n_shard) {
  fits[ss] <- density(model_out[[ss]][tslice,])
}


# consensus MCMC ensemble parameters



# calculate the weights for the individual componetns within each chain
weights <- vector('list',n_shard)
for (ss in 1:n_shard) {
  weights[[ss]] <- solve(cov(parameters[[ss]]))
}

# calculate the factor in big parentheses
summ <- NULL
for (ss in 1:n_shard) {
  summ <- sum + weights[[ss]]
}
winv <- solve(summ)

# combine at each iteration according to the weighting above
parametersC <- array(0, dim=dim(parameters[[1]]))
summ <- NULL
for (ss in 1:n_shard) {
  summ <- summ + weights[[ss]] %*% parameters[[ss]]
}

for (t in 1:nrow(parameters[[1]])) {
    ##parametersC[t,] <- winv %*% (weights1 %*% samples1[t,] + weights2 %*% samples2[t,])
    parametersC[t,] <- winv %*% (summ)
}

#===============================================================================
# End
#===============================================================================
