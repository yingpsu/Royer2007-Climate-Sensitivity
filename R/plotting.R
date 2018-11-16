##==============================================================================
## plotting.R
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')
plot.dir <- '../figures/'
load('../output/analysis.RData')

library(Hmisc)

##==============================================================================
# Figure .. Observations and fitted likelihood surface.

pdf(paste(plot.dir,'data_likelihood.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.75,.15,.15))
plot(-time, likelihood_quantiles[,'50'], type='l', xlim=c(-450,0), ylim=c(0,6500), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n')
polygon(-c(time,rev(time)), c(likelihood_quantiles[,1],rev(likelihood_quantiles[,9])), col='lightcyan1', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,2],rev(likelihood_quantiles[,8])), col='lightcyan2', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,3],rev(likelihood_quantiles[,7])), col='lightcyan3', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,4],rev(likelihood_quantiles[,6])), col='lightcyan4', border=NA)
lines(-time, likelihood_quantiles[,'50'], lwd=2)
points(-data_calib_ctrl$age, data_calib_ctrl$co2, pch='+', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.4, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=2.4, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()

##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

#model_quantiles[,'maxpost'] <- model_out[,which.max(lpost_out)]


## Log scale (model and points, no likelihood surface)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'q50']), type='l', xlim=c(-450,0), ylim=c(27,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='seagreen1', border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q05'],rev(model_quantiles[,'q95']))), col='gray', border=NA)
#lines(-time, log10(likelihood_quantiles[,'50']), lwd=2, lty=2)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
#lines(-time, model_ref, lwd=2, lty=2)
points(-data_calib_ctrl$age, log10(data_calib_ctrl$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-450, log10(50), c('Data','Max posterior','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, log10(50), c('Data','Max posterior','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


## Log scale (obs and likelihood surface, for SOM)
likelihood_nonansense <- cbind(time[which(!is.na(likelihood_quantiles[,2]))],
                               likelihood_quantiles[which(!is.na(likelihood_quantiles[,2])),2],
                               likelihood_quantiles[which(!is.na(likelihood_quantiles[,2])),8],
                               likelihood_quantiles[which(!is.na(likelihood_quantiles[,2])),'50'])
likelihood_nonansense[which(likelihood_nonansense[,2] <= 0),2] = 0.0001

pdf(paste(plot.dir,'likelihood_vs_obspts_logscale_SOM.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(likelihood_quantiles[,'50']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
#polygon(-c(time,rev(time)), log10(c(likelihood_quantiles[,1],rev(likelihood_quantiles[,9]))), col='lightcyan1', border=NA)
polygon(-c(likelihood_nonansense[,1],rev(likelihood_nonansense[,1])), log10(c(likelihood_nonansense[,2],rev(likelihood_nonansense[,3]))), col='gray', border=NA)
#polygon(-c(time,rev(time)), log10(c(likelihood_quantiles[,3],rev(likelihood_quantiles[,7]))), col='lightcyan3', border=NA)
#polygon(-c(time,rev(time)), log10(c(likelihood_quantiles[,4],rev(likelihood_quantiles[,6]))), col='lightcyan4', border=NA)
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='seagreen1', border=NA)
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q05'],rev(model_quantiles[,'q95']))), col='gray', border=NA)
lines(-likelihood_nonansense[,1], log10(likelihood_nonansense[,4]), lwd=2, lty=1, col='black')
#lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
#lines(-time, model_ref, lwd=2, lty=2)
points(-data_calib_ctrl$age, log10(data_calib_ctrl$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-450, log10(50), c('Data','Likelihood median','5-95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.9, bty='n')
legend(-450, log10(50), c('Data','Likelihood median','5-95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.9, bty='n')
dev.off()


##==============================================================================
# Figure 3. Posterior probability density for Earth system sensitivity parameter
# (deltaT2X), relative to previous studies.

#plot(deltaT2X_density$x, deltaT2X_density$y, type='l')

# get priors too
row_num <- match('deltaT2X',input$parameter)
x_cs <- seq(from=0, to=10, by=0.1)
f_cs <- dlnorm(x=x_cs, meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
row_num <- match('GLAC',input$parameter)
x_gl <- seq(from=0, to=10, by=0.1)
f_gl <- dnorm(x=x_gl, mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))

# Royer et al 2007:  1.5 and 6.2 deg C (5–95% range), 2.8 best fit
x_royer2007 <- c(1.6, 2.8, 5.5)
# Park and Royer 2011 from CSV/Excel table
x_pr2011 <- pr2011_cdf(c(.05,.5,.95))
iglac <- match('GLAC',parnames_calib)
x_thisstudy <- quantile(parameters[,ics], c(.05,.5,.95))  # 3.148759 4.203719 5.621429
x_glac <- quantile(parameters[,ics]*parameters[,iglac], c(.05,.5,.95))  # 5.678216  9.289884 14.268706
x_norm <- quantile(parameters_nm[,ics], c(.05,.5,.95))  # 3.273916 4.396285 5.846267
x_ktc2017 <- c(3.7, 5.6, 7.5)

offset <- 0.06

pdf(paste(plot.dir,'deltaT2X.pdf',sep=''),width=4,height=3, colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=2, xlim=c(0.9,10.5), ylim=c(0,.7+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
#polygon(-c(time,rev(time)), c(model_quantiles[,'q025'],rev(model_quantiles[,'q975'])), col='aquamarine1', border=NA)
#polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='aquamarine3', border=NA)
#lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=2, lty=3)
lines(deltaT2X_density_pr2011$x, deltaT2X_density_pr2011$y + offset, lwd=2, lty=2)
lines(x_cs, f_cs + offset, lwd=2, lty=3)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.2, cex=1)
mtext('Density', side=2, line=0.3, cex=1)
arrows(1, 0, 1, .7+offset, length=0.08, angle=30, code=2)
axis(1, at=seq(0,10), cex.axis=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)
y0 <- 0.7*offset; arrows(x_thisstudy[1], y0, x_thisstudy[3], y0, length=0.05, angle=90, code=3); points(x_thisstudy[2], y0, pch=16)
#y1 <- 0.35*offset; arrows(x_royer2007[1], y1, x_royer2007[3], y1, length=0.05, angle=90, code=3); points(x_royer2007[2], y1, pch=15)
y1 <- 0.35*offset; arrows(x_pr2011[1], y1, x_pr2011[3], y1, length=0.05, angle=90, code=3); points(x_pr2011[2], y1, pch=15)
#y2 <- 0.08; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.05, angle=90, code=3); points(x_ktc2017[2], y2, pch=17)
legend(4.81,0.8, c('5-95% range, PR2011','5-95% range, this study','Posterior, PR2011','Posterior, this study','Prior, both studies'),
       pch=c(15,16,NA,NA,NA), lty=c(1,1,2,1,3), cex=.95, bty='n')
dev.off()


offset <- 0.08

pdf(paste(plot.dir,'deltaT2X_SOM.pdf',sep=''),width=6,height=4, colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=2, xlim=c(0.8,20), ylim=c(0,.63+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=2, lty=2)
#lines(deltaT2X_density_pr2011$x, deltaT2X_density_pr2011$y + offset, lwd=2, lty=2)
lines(c(10,20), c(offset,offset), lty=1, lwd=2)
#lines(x_cs, f_cs + offset, lwd=2, lty=3)
lines(deltaT2Xglac_density$x, offset+deltaT2Xglac_density$y, lwd=2, lty=3)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.2, cex=1)
mtext('Density', side=2, line=0.3, cex=1)
arrows(1, 0, 1, .6+offset, length=0.08, angle=30, code=2)
axis(1, at=seq(0,20,1), labels=rep('',21), col='gray')
axis(1, at=seq(0,20,5), labels=c('0','5','10','15','20'), cex.axis=1)
y0 <- 0.16*offset; arrows(x_thisstudy[1], y0, x_thisstudy[3], y0, length=0.05, angle=90, code=3); points(x_thisstudy[2], y0, pch=16)
y1 <- 0.39*offset; arrows(x_norm[1], y1, x_norm[3], y1, length=0.05, angle=90, code=3); points(x_norm[2], y1, pch=2)
y3 <- 0.85*offset; arrows(x_glac[1], y3, x_glac[3], y3, length=0.05, angle=90, code=3); points(x_glac[2], y3, pch=1)
y2 <- 0.62*offset; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.05, angle=90, code=3); points(x_ktc2017[2], y2, pch=17)
legend(7,0.7, c('5-95% range, KTC2017','5-95% range, this study','5-95% range (normal assumption), this study','5-95% range (glacial), this study',
                'Posterior, this study','Posterior (normal assumption), this study','Posterior (glacial), this study'), pch=c(17,16,2,1,NA,NA,NA), lty=c(1,1,1,1,2,3), cex=.95, bty='n')
dev.off()



print(paste('deltaT2X 5-50-95% quantiles =',quantile(parameters[,ics], c(.05,.5,.95))))
print(paste('deltaT2X 16-50-84% quantiles =',quantile(parameters[,ics], c(.16,.5,.84))))
print(paste('... using normally distributed errors =',quantile(parameters_nm[,ics], c(.05,.5,.95))))
print(paste('fraction of dT2X >= 6 is:',length(which(parameters[,ics]>=6))/nrow(parameters)))
print(paste('GLAC 5-50-95% quantiles =',quantile(parameters[,iglac], c(.05,.5,.95))))
print(paste('glacial dT2X 5-50-95% quantiles =',quantile(parameters[,iglac]*parameters[,ics], c(.05,.5,.95))))
print(paste('glacial dT2X 16-50-84% quantiles =',quantile(parameters[,iglac]*parameters[,ics], c(.16,.5,.84))))
print(paste('PR2011 dT2X 5-50-95% quantiles =',pr2011_cdf(c(.05,.5,.95))))
print(paste('PR2011 dT2X 16-50-84% quantiles =',pr2011_cdf(c(.16,.5,.84))))

##==============================================================================
# Figure 4. Radial convergence diagrams for the sensitivity experiment. The
# “Sensitive” parameters are those identified as significant in the sub-ensemble
# correlation evaluation. Filled blue nodes represent first-order sensitivity
# indices; filled purple nodes represent total-order sensitivity indices; filled
# gray bars represent second-order sensitivity indices for the interaction
# between the parameter pair.

# TODO!!



##==============================================================================
# Figure S2.  Evidence of multi-modality, and dispersion of probability for
# higher CO2 data points

pdf(paste(plot.dir,'likelihoodslice_SOM.pdf',sep=''),width=4,height=3, colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.7,.3,.1,.25))
plot(mm_example$co2, mm_example$fit, type='l', lwd=2, xlim=c(-70,5000), ylim=c(0,6.5e-4),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
text(3000, 5.5e-4, paste('Likelihood surface slice\n at age =',mm_example$age,'Mya'))
axis(1, at=seq(0,5000,200), labels=rep('',length(seq(0,5000,200))), col='gray')
axis(1, at=seq(0,5000,1000), labels=c('0','1000','2000','3000','4000','5000'), cex.axis=1)
arrows(0, 0, 0, 6.3e-4, length=0.08, angle=30, code=2)
mtext(expression(CO[2]~(ppmv)), side=1, line=2.4, cex=1)
mtext('Density', side=2, line=0.3, cex=1)
dev.off()


# check model ensemble during this time slice
fit <- density(model_out[34,])
mm_modeled <- mm_example
mm_modeled$fit <- fit$y
mm_modeled$co2 <- fit$x
plot(mm_modeled$co2, mm_modeled$fit, type='l', lwd=2, lty=2, xlim=c(0,5000),
     xlab='CO2 (ppmv)', ylab='Probability density')
lines(mm_example$co2, mm_example$fit, lwd=2)


# check precalibration ensemble during this time slice
parameters_precal <- readRDS('../output/precal_parameters_N2M-BS10k_13Nov2018.rds')
model_out_precal <- sapply(X=1:nrow(parameters_precal),
              FUN=function(k){model_forMCMC(par_calib=parameters_precal[k,],
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
fit2 <- density(model_out_precal[34,])
mm_precal <- mm_example
mm_precal$co2 <- fit2$x
mm_precal$fit <- fit2$y


# check parameters straight from the priors
library(lhs)
## draw parameters by Latin Hypercube (Sample)
parameters_lhs0 <- randomLHS(100000, n_parameters)

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

model_out_priors <- sapply(1:nrow(parameters_lhs), function(ss) {
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
fit3 <- density(model_out_priors[34,which(model_out_priors[34,] < 1e4)])
mm_priors <- mm_example
mm_priors$co2 <- fit3$x
mm_priors$fit <- fit3$y


plot(mm_modeled$co2, mm_modeled$fit, type='l', lwd=2, lty=2, xlim=c(0,5000),
     xlab='CO2 (ppmv)', ylab='Probability density')
lines(mm_example$co2, mm_example$fit, lwd=2)
lines(mm_precal$co2, mm_precal$fit, lwd=2, lty=3)
lines(mm_priors$co2, mm_priors$fit, lwd=2, lty=4)


##==============================================================================
# Figure S3. Posterior probability density for Earth system sensitivity
# parameter (deltaT2X), relative to previous studies), assuming a symmetric
# (Gaussian) error structure for the proxy data as opposed to skew-normal.

deltaT2X_density_nm <- density(parameters_nm[,ics], from=0, to=10)

#plot(deltaT2X_density$x, deltaT2X_density$y, type='l', lwd=2); lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y, type='l', lty=2, lwd=2)



##==============================================================================
# Other scratch work

if(FALSE) {


par(mfrow=c(7,8))
for (p in 1:length(parnames)) {plot(parameters[,p], type='l', ylab=parnames[p])}

par(mfrow=c(1,1))
hist(parameters[round(nrow(parameters)*0.5):nrow(parameters), match('deltaT2X',parnames)],
     freq=FALSE, xlab='deltaT2X [deg C]', main='')


}
# End other scratch work

##==============================================================================


##==============================================================================
## End
##==============================================================================
