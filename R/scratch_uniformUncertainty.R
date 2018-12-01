#===============================================================================
#
#===============================================================================

setwd('~/codes/GEOCARB/output')
load('GEOCARB_MCMC_unifUnc_28Nov2018nm-unifUnc.RData')
library(Hmisc)

# histograms figure
par(mfrow=c(5,4), mai=c(.3,.35,.1,.01))
for (mm in 1:length(amcmc_par1)) {hist(amcmc_par1[[mm]]$samples[5e5:1e6,10], freq=FALSE, xlim=c(1,9), main='', breaks=seq(0,10,by=0.25))}

# history plots figure
par(mfrow=c(5,4), mai=c(.3,.35,.1,.01))
for (mm in 1:length(amcmc_par1)) {plot(amcmc_par1[[mm]]$samples[,10], type='l')}

# grab posterior parameter ensembles from each shard
parameters <- amcmc_par[[1]]$samples[seq((niter_mcmc*0.6+1),niter_mcmc,by=100),]
for (mm in 2:length(amcmc_par1)) {
    parameters <- rbind(parameters, amcmc_par1[[mm]]$samples[seq((niter_mcmc*0.6+1),niter_mcmc,by=100),])
}

# run posterior model ensembles for each shard
source('../R/run_geocarbF.R')
model_out <- sapply(1:nrow(parameters), function(ii) {
    model_forMCMC(par_calib=parameters[ii,],
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

# posterior time slice at t=34 (age=240 Myr) from each shard
tslice <- 34
fits <- density(model_out[tslice,])


# get quantiles
model_quantiles <- mat.or.vec(nr=n_time, nc=3)
colnames(model_quantiles) <- c('q05','q50','q95')
for (t in 1:n_time) {
  model_quantiles[t,1:3] <- quantile(model_out[t,], c(.05,.50,.95))
}


# plot of 240 Myr time slice from each shard, and consensus

par(mfrow=c(2,1), mai=c(.9,.9,.1,.1))
plot(fits$x, fits$y, type='l', lty=2, lwd=2,
     xlab='CO2 (ppmv)', ylab='Probability density', xlim=c(0,5000), ylim=c(0,0.0085))
lines(mm_example$co2, mm_example$fit, lty=1, lwd=2)
legend(1000,0.008, c('3 separate MCMC chains','Full likelihood surface'), lty=c(2,1), lwd=2, bty='n')



# plot of hindcast CO2
pdf(paste('../figures/model_ensemble_vs_obspts_unifUnc.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.6,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q05'],rev(model_quantiles[,'q95']))), col='gray', border=NA)
lines(-time, log10(model_quantiles[,'q50']), lwd=2)
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

plot(-time, model_quantiles[,'q50'], type='l', xlim=c(-450,0), ylim=c(0,(4100)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='gray', border=NA)
lines(-time, model_quantiles[,'q50'], lwd=2)
points(-data_calib$age, data_calib$co2, pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
axis(2, at=seq(0,4000,500), labels=seq(0,4000,500), cex.axis=1.1, las=1)
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

dev.off()


#===============================================================================

# comparison against precalibration parameters
parameters_precal <- readRDS('precal_parameters_N2M-BS10k_13Nov2018.rds')

# look at the likelihood function, priors, posterior for the simulations that
# hit the high peak around 200 Myr and those that hit the low peak

model_precal <- sapply(1:nrow(parameters_precal), function(ii) {
    model_forMCMC(par_calib=parameters_precal[ii,],
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

model_quantiles_precal <- mat.or.vec(nr=n_time, nc=3)
colnames(model_quantiles_precal) <- c('q05','q50','q95')
for (t in 1:n_time) {
  model_quantiles_precal[t,1:3] <- quantile(model_precal[t,], c(.05,.50,.95))
}

pdf(paste('../figures/model_ensemble_vs_obspts_precal.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles_precal[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.6,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles_precal[,'q05'],rev(model_quantiles_precal[,'q95']))), col='gray', border=NA)
lines(-time, log10(model_quantiles_precal[,'q50']), lwd=2)
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

plot(-time, model_quantiles_precal[,'q50'], type='l', xlim=c(-450,0), ylim=c(0,(4100)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), c(model_quantiles_precal[,'q05'],rev(model_quantiles_precal[,'q95'])), col='gray', border=NA)
lines(-time, model_quantiles_precal[,'q50'], lwd=2)
points(-data_calib$age, data_calib$co2, pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
axis(2, at=seq(0,4000,500), labels=seq(0,4000,500), cex.axis=1.1, las=1)
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

dev.off()

#===============================================================================

#
# check parameters straight from the priors
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

model_priors <- sapply(1:nrow(parameters_lhs), function(ss) {
                    model_forMCMC(par_calib=parameters_lhs[ss,],
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

model_quantiles_priors <- mat.or.vec(nr=n_time, nc=3)
colnames(model_quantiles_priors) <- c('q05','q50','q95')
for (t in 1:n_time) {
  model_quantiles_priors[t,1:3] <- quantile(model_priors[t,], c(.05,.50,.95))
}

pdf(paste('../figures/model_ensemble_vs_obspts_priors.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles_priors[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.6,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
#polygon(-c(time,rev(time)), log10(c(model_quantiles_priors[,'q05'],rev(log10(model_quantiles_priors[,'q95'])))), col='gray', border=NA)
polygon(c(-time[-58],rev(-time)[-1]), c(log10(model_quantiles_priors[-58,'q95']),rev(log10(model_quantiles_priors[,'q05']))[-1]), col='gray', border=NA)
lines(-time, log10(model_quantiles_priors[,'q50']), lwd=2)
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

plot(-time, model_quantiles_priors[,'q50'], type='l', xlim=c(-450,0), ylim=c(0,(4100)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(c(-time[-58],rev(-time)[-1]), c(model_quantiles_priors[-58,'q95'],rev(model_quantiles_priors[,'q05'])[-1]), col='gray', border=NA)
lines(-time, model_quantiles_priors[,'q50'], lwd=2)
points(-data_calib$age, data_calib$co2, pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
axis(2, at=seq(0,4000,500), labels=seq(0,4000,500), cex.axis=1.1, las=1)
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, 50, c('Data','Median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

dev.off()

#===============================================================================

# check what the set of simulations that hit the upper mode looks like, and
# what the set that hits the lower mode looks like

model_ensemble <- model_precal

idx_lo <- which(model_ensemble[34,] <= 1610)
idx_hi <- which(model_ensemble[34,] > 1610)

model_quantiles_hi <- model_quantiles_lo <- mat.or.vec(nr=n_time, nc=3)
colnames(model_quantiles_hi) <- colnames(model_quantiles_lo) <- c('q05','q50','q95')
for (t in 1:n_time) {
  model_quantiles_hi[t,1:3] <- quantile(model_ensemble[t,idx_hi], c(.05,.50,.95))
  model_quantiles_lo[t,1:3] <- quantile(model_ensemble[t,idx_lo], c(.05,.50,.95))
}

pdf(paste('../figures/model_ensemble_vs_obspts_hilo.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles_hi[,'q50']), type='l', xlim=c(-450,0), ylim=c(0.6,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles_hi[,'q05'],rev(model_quantiles_hi[,'q95']))), col=rgb(1,0,0,.3), border=NA)
lines(-time, log10(model_quantiles_hi[,'q50']), lwd=2, col=rgb(1,0,0))
polygon(-c(time,rev(time)), log10(c(model_quantiles_lo[,'q05'],rev(model_quantiles_lo[,'q95']))), col=rgb(0,1,0,.3), border=NA)
lines(-time, log10(model_quantiles_lo[,'q50']), lwd=2, col=rgb(0,1,0))
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

dev.off()


#===============================================================================


# check how the simulations that hit the top compare

lpri_con <- sapply(1:nrow(parametersC), function(ss) {
             log_prior(par_calib=parametersC[ss,],
                        par_fixed=par_fixed0,
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
                        bounds_calib=bounds_calib,
                        do_sample_tvq=DO_SAMPLE_TVQ)})

llike_con <- sapply(1:nrow(parametersC), function(ss) {
              log_like(par_calib=parametersC[ss,],
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
                        likelihood_fit=likelihood_fit,
                        idx_data=idx_data,
                        do_sample_tvq=DO_SAMPLE_TVQ,
                        par_time_center=par_time_center,
                        par_time_stdev=par_time_stdev)})

idx_hic <- (which(model_con[34,] > 1610))
idx_loc <- (which(model_con[34,] <= 1610))

like_hi <- density(llike_con[idx_hic])
like_lo <- density(llike_con[idx_loc])
pri_hi <- density(lpri_con[idx_hic])
pri_lo <- density(lpri_con[idx_loc])

par(mfrow=c(2,1))
plot(like_lo$x, like_lo$y, xlim=c(-320,-280), type='l', lty=1, xlab='log-probability', ylab='density')
  lines(like_hi$x, like_hi$y, lty=2)
  legend(-320,0.15,c('Low mode','High mode'), lty=c(1,2), bty='n')
plot(pri_lo$x, pri_lo$y, xlim=c(-33, 2), ylim=c(0,0.12), type='l', xlab='log-probability', ylab='density')
  lines(pri_hi$x, pri_hi$y)

#===============================================================================
# End
#===============================================================================
