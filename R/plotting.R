##==============================================================================
## plotting.R
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')
plot.dir <- '../figures/'
load('../output/analysis.RData')

##==============================================================================
# Figure 1. Observations and fitted likelihood surface.

pdf(paste(plot.dir,'data_likelihood.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.75,.15,.15))
plot(-time, likelihood_quantiles[,'50'], type='l', xlim=c(-450,0), ylim=c(0,6500), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n')
polygon(-c(time,rev(time)), c(likelihood_quantiles[,1],rev(likelihood_quantiles[,9])), col='lightcyan1', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,2],rev(likelihood_quantiles[,8])), col='lightcyan2', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,3],rev(likelihood_quantiles[,7])), col='lightcyan3', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,4],rev(likelihood_quantiles[,6])), col='lightcyan4', border=NA)
lines(-time, likelihood_quantiles[,'50'], lwd=2)
points(-data_calib$age, data_calib$co2, pch='+', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.4, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=2.4, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
dev.off()

# TODO -- need to sort out why the early data are not plotting full polygons
# TODO -- need to sort out why the early data are not plotting full polygons
# TODO -- need to sort out why the early data are not plotting full polygons


##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

#model_quantiles[,'maxpost'] <- model_out[,which.max(lpost_out)]


# TODO -- switch from plotting against hte obs directly to plotting against the fitted likelihood surface (above)
# TODO -- switch from plotting against hte obs directly to plotting against the fitted likelihood surface (above)
# TODO -- switch from plotting against hte obs directly to plotting against the fitted likelihood surface (above)


pdf(paste(plot.dir,'model_ensemble_vs_obs.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.75,.15,.15))
plot(-time, model_quantiles[,'q50'], type='l', xlim=c(-450,0), ylim=c(0,6500), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n')
polygon(-c(time,rev(time)), c(likelihood_quantiles[,1],rev(likelihood_quantiles[,9])), col='lightcyan1', border=NA)
polygon(-c(time,rev(time)), c(likelihood_quantiles[,2],rev(likelihood_quantiles[,8])), col='lightcyan2', border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[,'q025'],rev(model_quantiles[,'q975'])), col='seagreen1', border=NA)
polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='seagreen3', border=NA)
lines(-time, likelihood_quantiles[,'50'], lwd=2, lty=2)
lines(-time, model_quantiles[,'maxpost'], lwd=2)
#lines(-time, model_ref, lwd=2, lty=2)
#points(-data_calib$age, data_calib$co2, pch='+', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.4, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=2.4, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
dev.off()


##==============================================================================
# Figure 3. Posterior probability density for Earth system sensitivity parameter
# (deltaT2X), relative to previous studies.

#plot(deltaT2X_density$x, deltaT2X_density$y, type='l')

# Royer et al 2007:  1.5 and 6.2 deg C (5–95% range), 2.8 best fit
x_royer2007 <- c(1.5, 2.8, 6.2)
x_ktc2017 <- c(3.7, 5.6, 7.5)

offset <- 0.15

pdf(paste(plot.dir,'deltaT2X.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.8,.45,.15,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=2, xlim=c(1,9), ylim=c(0,.65+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
#polygon(-c(time,rev(time)), c(model_quantiles[,'q025'],rev(model_quantiles[,'q975'])), col='aquamarine1', border=NA)
#polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='aquamarine3', border=NA)
lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=2, lty=2)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.4, cex=1)
mtext('Density', side=2, line=0.85, cex=1)
axis(1, at=seq(0,10,0.5), labels=c('0','','1','','2','','3','','4','','5','','6','','7','','8','','9','','10'), cex.axis=1.1)
y1 <- 0.12; arrows(x_royer2007[1], y1, x_royer2007[3], y1, length=0.05, angle=90, code=3); points(x_royer2007[2], y1, pch=16)
y2 <- 0.08; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.05, angle=90, code=3); points(x_ktc2017[2], y2, pch=15)
legend(6,0.8, c('R2007','KTC2017'), pch=c(16,15), cex=1, bty='n')
dev.off()


##==============================================================================
# Figure 4. Radial convergence diagrams for the sensitivity experiment. The
# “Sensitive” parameters are those identified as significant in the sub-ensemble
# correlation evaluation. Filled blue nodes represent first-order sensitivity
# indices; filled purple nodes represent total-order sensitivity indices; filled
# gray bars represent second-order sensitivity indices for the interaction
# between the parameter pair.

# TODO!!



##==============================================================================
# Figure S2.  Posterior model ensemble (gray shaded region denotes 5-95%
# credible range), maximum posterior score simulation (solid bold line) and
# uncalibrated model simulation (dashed line), with proxy data points
# superimposed (+ markers), assuming a symmetric (Gaussian) error structure for
# the proxy data as opposed to skew-normal (main text).

#model_quantiles_nm



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
