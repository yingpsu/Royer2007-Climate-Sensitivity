#===============================================================================
#
#===============================================================================

setwd('~/codes/GEOCARB/output')
#load('GEOCARB_MCMC_tvq_split_16Nov2018sn-split.RData')
load('GEOCARB_MCMC_tvq_split_01Dec2018sn-split.RData')
library(Hmisc)

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
  fits[[ss]] <- density(model_out[[ss]][tslice,])
}


# consensus MCMC ensemble parameters


# CHECK THESE CALCULATIONS
# -- probably not working quite right because the
#    distribution for deltaT2X shifts farther to the left than any one of the
#    constituent distributions
# CHECK THESE CALCULATIONS



# calculate the weights for the individual componetns within each chain
weights <- vector('list',n_shard)
for (ss in 1:n_shard) {
  weights[[ss]] <- solve(cov(parameters[[ss]]))
}

# calculate the factor in big parentheses
summ <- weights[[1]]
for (ss in 2:n_shard) {
  summ <- summ + weights[[ss]]
}
winv <- solve(summ)

# combine at each iteration according to the weighting above
summ <- weights[[1]] %*% t(parameters[[1]])
for (ss in 2:n_shard) {
  summ <- summ + weights[[ss]] %*% t(parameters[[ss]])
}

parametersC <- array(0, dim=dim(parameters[[1]]))
parametersC <- t(winv %*% (summ))
#for (t in 1:nrow(parameters[[1]])) {
#    parametersC[t,] <- winv %*% (weights1 %*% samples1[t,] + weights2 %*% samples2[t,])
#}

# run posterior model ensemble using consensus samples
model_con <- sapply(1:nrow(parametersC), function(ii) {
        model_forMCMC(par_calib=parametersC[ii,],
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

fitsC <- density(model_con[tslice,])


# get quantiles
model_quantilesC <- mat.or.vec(nr=n_time, nc=3)
colnames(model_quantilesC) <- c('q05','q50','q95')
for (t in 1:n_time) {
  idx <- which(!is.infinite(model_con[t,]))
  model_quantilesC[t,1:3] <- quantile(model_con[t,idx], c(.05,.50,.95), na.rm=TRUE)
}


# plot of 240 Myr time slice from each shard, and consensus

par(mfrow=c(2,1), mai=c(.9,.9,.1,.1))
plot(fits[[1]]$x, fits[[1]]$y, type='l', lty=2, lwd=2,
     xlab='CO2 (ppmv)', ylab='Probability density', xlim=c(0,5000), ylim=c(0,0.0085))
lines(fits[[2]]$x, fits[[2]]$y, lty=2, lwd=2)
lines(fits[[3]]$x, fits[[3]]$y, lty=2, lwd=2)
lines(mm_example$co2, mm_example$fit, lty=1, lwd=2)
legend(1000,0.008, c('3 separate MCMC chains','Full likelihood surface'), lty=c(2,1), lwd=2, bty='n')


# plot of distributions of ESS from each shard, and consensus

fitp <- vector('list',n_shard)
for (ss in 1:n_shard) {fitp[[ss]] <- density(parameters[[ss]][,10])}
fitp_con <- density(parametersC[,10])

plot(fitp[[1]]$x, fitp[[1]]$y, type='l', lty=2, lwd=2,
     xlab='deltaT2X (deg C)', ylab='Probability density', xlim=c(1,10))
lines(fitp[[2]]$x, fitp[[2]]$y, lty=2, lwd=2)
lines(fitp[[3]]$x, fitp[[3]]$y, lty=2, lwd=2)
#lines(mm_example$co2, mm_example$fit, lty=1, lwd=2)
legend(1000,0.008, c('3 separate MCMC chains','Full likelihood surface'), lty=c(2,1), lwd=2, bty='n')



# plot of hindcast CO2
pdf(paste('../figures/model_ensemble_vs_obspts_consensus.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantilesC[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.6,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantilesC[,'q05'],rev(model_quantilesC[,'q95']))), col='gray', border=NA)
lines(-time, log10(model_quantilesC[,'q50']), lwd=2)
points(-data_calib$age, log10(data_calib$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-450, log10(50), c('Data','Median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, log10(50), c('Data','Median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

plot(-time, model_quantilesC[,'q50'], type='l', xlim=c(-450,0), ylim=c(0,(4100)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), c(model_quantilesC[,'q05'],rev(model_quantilesC[,'q95'])), col='gray', border=NA)
lines(-time, model_quantilesC[,'q50'], lwd=2)
points(-data_calib$age, data_calib$co2, pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
axis(2, at=seq(0,4000,500), labels=seq(0,4000,500), cex.axis=1.1, las=1)
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

dev.off()

#####
#####
#####

library(lhs)
parameters_lhs0 <- randomLHS(100000, 68)

## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
parameters_lhs <- parameters_lhs0  # initialize
colnames(parameters_lhs) <- parnames_calib
for (i in 1:n_const_calib) {
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    parameters_lhs[,i] <- qnorm(p=parameters_lhs0[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    parameters_lhs[,i] <- qlnorm(p=parameters_lhs0[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}
for (i in (n_const_calib+1):length(parnames_calib)) {
  parameters_lhs[,i] <- qbeta(p=parameters_lhs0[,i], shape1=5, shape2=5)
}

ss <- 3
tmp <- sapply(1:nrow(parameters_lhs), function(k) {
               log_like(par_calib=parameters_lhs[k,],
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
                        data_calib=data_calib,
                        ind_mod2obs=ind_mod2obs,
                        ind_expected_time=ind_expected_time,
                        ind_expected_const=ind_expected_const,
                        iteration_threshold=iteration_threshold,
                        loglikelihood_smoothed=loglikelihood_smoothed,
                        likelihood_fit=likelihood_fit[[ss]],
                        idx_data=idx_data[[ss]],
                        do_sample_tvq=DO_SAMPLE_TVQ,
                        par_time_center=par_time_center,
                        par_time_stdev=par_time_stdev)})


#===============================================================================
# End
#===============================================================================
