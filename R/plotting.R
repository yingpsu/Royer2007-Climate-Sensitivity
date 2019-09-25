##==============================================================================
## plotting.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')
plot.dir <- '../figures/'
load('../output/analysis_24Sep2019.RData')

library(Hmisc)

##==============================================================================
# Figure 5 (Methods). Observations and fitted likelihood surface.

source('likelihood_surface_quantiles.R')

pdf(paste(plot.dir,'data_likelihood.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(likelihood_quantiles[,'50']), type='l', xlim=c(-430,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'025'],rev(likelihood_quantiles[idx_likelihood_ages,'975']))), col='gray85', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'05'],rev(likelihood_quantiles[idx_likelihood_ages,'95']))), col='gray70', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'17'],rev(likelihood_quantiles[idx_likelihood_ages,'83']))), col='gray55', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'25'],rev(likelihood_quantiles[idx_likelihood_ages,'75']))), col='gray40', border=NA)
lines(-time[idx_likelihood_ages], log10(likelihood_quantiles[idx_likelihood_ages,'50']), lwd=2, lty=1)
points(-data_calib$age, log10(data_calib$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-452, log10(40), c('Data','Likelihood median','95% range','90% range','66% range','50% range'), pch=c(4,NA,15,15,15,15), col=c('black','black','gray85','gray70','gray55','gray40'), cex=.8, bty='n')
legend(-452, log10(40), c('Data','Likelihood median','95% range','90% range','66% range','50% range'), pch=c(NA,'-',NA,NA,NA,NA), col=c('black','black','gray85','gray70','gray55','gray40'), cex=.8, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()

# version with the maximum posterior score model simulation superimposed too
pdf(paste(plot.dir,'data_likelihood_withModel.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(likelihood_quantiles[,'50']), type='l', xlim=c(-430,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'025'],rev(likelihood_quantiles[idx_likelihood_ages,'975']))), col='gray85', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'05'],rev(likelihood_quantiles[idx_likelihood_ages,'95']))), col='gray70', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'17'],rev(likelihood_quantiles[idx_likelihood_ages,'83']))), col='gray55', border=NA)
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'25'],rev(likelihood_quantiles[idx_likelihood_ages,'75']))), col='gray40', border=NA)
lines(-time[idx_likelihood_ages], log10(likelihood_quantiles[idx_likelihood_ages,'50']), lwd=2, lty=1)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2, lty=5, col="salmon3")
points(-data_calib$age, log10(data_calib$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-452, log10(40), c('Data','Likelihood median','95% range','90% range','66% range','50% range'), pch=c(4,NA,15,15,15,15), col=c('black','black','gray85','gray70','gray55','gray40'), cex=.8, bty='n')
legend(-452, log10(40), c('Data','Likelihood median','95% range','90% range','66% range','50% range'), pch=c(NA,'-',NA,NA,NA,NA), col=c('black','black','gray85','gray70','gray55','gray40'), cex=.8, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

#model_quantiles[,'maxpost'] <- model_out[,which.max(lpost_out)]
idx_gastaldo <- which(data_calib_all$reference=="Gastaldo et al., 2014")


## Log scale (model and points, no likelihood surface)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_errbars.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='gray', border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='orange', pch=16, cex=1)
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='black', pch=16, cex=0.5)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.75, bty='n')
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


## Log scale (model and points, with CO2 error bars, no likelihood surface)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_errbars.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='gray', border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='orange', pch=16, cex=1)
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='black', pch=16, cex=0.5)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.75, bty='n')
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


## Log scale (model and points, with CO2 and age error bars, no likelihood surface)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_2errbars.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='gray', border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='orange', pch=16, cex=1)
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='black', pch=16, cex=0.5)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.75, bty='n')
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


## Log scale (model and points, with CO2 and age error bars, no likelihood surface, but with PR2011 too)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_2errbars.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='gray', border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='orange', pch=16, cex=1)
points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='black', pch=16, cex=0.5)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.75, bty='n')
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


pdf(paste(plot.dir,'model_ensemble_vs_royer_logscale.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
points(-data_calib$age, log10(data_calib$co2), pch=4, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
polygon(-c(time,rev(time)), log10(c(model_quantiles_royer[,'q025'],rev(model_quantiles_royer[,'q975']))), col=rgb(.6,.2,.6,.5), border=NA)
lines(-time, log10(model_quantiles_royer[,'co2']), lwd=2, lty=5, col=rgb(.6,.2,.6))
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.2,.6,.6,.5), border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2, col=rgb(.2,.6,.6))
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(4,NA,15), col=c('black','black','gray'), cex=.75, bty='n')
legend(-457, log10(35), c('Data','Max posterior','95% range'), pch=c(NA,'-',NA), col=c('black','black','gray'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()

##==============================================================================



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

pdf(paste(plot.dir,'deltaT2X.pdf',sep=''),width=4,height=3, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=1.7, xlim=c(0.9,10.1), ylim=c(0,.9+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE, col="steelblue")
#polygon(-c(time,rev(time)), c(model_quantiles[,'q025'],rev(model_quantiles[,'q975'])), col='aquamarine1', border=NA)
#polygon(-c(time,rev(time)), c(model_quantiles[,'q05'],rev(model_quantiles[,'q95'])), col='aquamarine3', border=NA)
#lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=2, lty=3)
lines(deltaT2X_density_pr2011$x, deltaT2X_density_pr2011$y + offset, lwd=1.7, lty=2)
lines(x_cs, f_cs + offset, lwd=1.7, lty=3, col="steelblue")
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.3)
mtext('Density', side=2, line=0.3)
arrows(1, 0, 1, .85+offset, length=0.08, angle=30, code=2)
axis(1, at=seq(0,10))
minor.tick(nx=4, ny=0, tick.ratio=0.5)
y0 <- 0.7*offset; arrows(x_thisstudy[1], y0, x_thisstudy[3], y0, lwd=1.7, length=0.05, angle=90, code=3, col="steelblue"); points(x_thisstudy[2], y0, pch=16, col="steelblue")
#y1 <- 0.35*offset; arrows(x_royer2007[1], y1, x_royer2007[3], y1, length=0.05, angle=90, code=3); points(x_royer2007[2], y1, pch=15)
y1 <- 0.35*offset; arrows(x_pr2011[1], y1, x_pr2011[3], y1, length=0.05, angle=90, code=3); points(x_pr2011[2], y1, pch=15)
#y2 <- 0.08; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.05, angle=90, code=3); points(x_ktc2017[2], y2, pch=17)
legend(4.93,1, c('5-95% range, PR2011','PR2011','5-95% range, this study','Posterior, this study','Prior, both studies'),
       pch=c(15,NA,16,NA,NA), lty=c(1,2,1,1,3), col=c("black","black","steelblue","steelblue","steelblue"), bty='n', lwd=1.7, cex=0.9)
dev.off()


##==============================================================================
# Figure S3. Posterior probability density for Earth system sensitivity
# parameter (deltaT2X), relative to previous studies), assuming a symmetric
# (Gaussian) error structure for the proxy data as opposed to skew-normal.

offset <- 0.08

pdf(paste(plot.dir,'deltaT2X_SOM.pdf',sep=''),width=6,height=4, colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.3,.13,.15))
plot(deltaT2X_density$x, deltaT2X_density$y + offset, type='l', lwd=1.7, xlim=c(0.8,20), ylim=c(0,.63+offset),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE, col='steelblue')
lines(deltaT2X_density_nm$x, deltaT2X_density_nm$y + offset, lwd=1.7, lty=2, col='seagreen3')
lines(c(10,20), c(offset,offset), lty=1, lwd=1.7)
lines(deltaT2Xglac_density$x, offset+deltaT2Xglac_density$y, lwd=1.7, lty=4, col='salmon3')
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.4)
mtext('Density', side=2, line=0.3)
arrows(1, 0, 1, .6+offset, length=0.08, angle=30, code=2)
axis(1, at=seq(0,20,1), labels=rep('',21), col='gray')
axis(1, at=seq(0,20,5), labels=c('0','5','10','15','20'), cex.axis=1)
y0 <- 0.16*offset; arrows(x_thisstudy[1], y0, x_thisstudy[3], y0, length=0.05, angle=90, code=3, lwd=1.5, col='steelblue'); points(x_thisstudy[2], y0, pch=16, col='steelblue')
y1 <- 0.39*offset; arrows(x_norm[1], y1, x_norm[3], y1, length=0.05, angle=90, code=3, lwd=1.5, lty=5, col='seagreen3'); points(x_norm[2], y1, pch=2, col='seagreen3')
y3 <- 0.85*offset; arrows(x_glac[1], y3, x_glac[3], y3, length=0.05, angle=90, code=3, lwd=1.5, lty=4, col='salmon3'); points(x_glac[2], y3, pch=1, col='salmon3')
y2 <- 0.62*offset; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.05, angle=90, lwd=1.5, code=3); points(x_ktc2017[2], y2, pch=17)
legend(7,0.7, c('5-95% range, Krissansen-Totton & Catling (2017)','5-95% range, this study','5-95% range (normal assumption), this study',
                '5-95% range (glacial), this study', 'Posterior, this study','Posterior (normal assumption), this study',
                'Posterior (glacial), this study'), pch=c(17,16,2,1,NA,NA,NA), lty=c(1,1,5,3,1,2,4), cex=.9, bty='n', lwd=1.5,
       col=c('black','steelblue','seagreen3','salmon3','steelblue','seagreen3','salmon3'))
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

# run plotting_sobol.R
source('plotting_sobol.R')


##==============================================================================
# Figure S2.  Evidence of multi-modality, and dispersion of probability for
# higher CO2 data points

pdf(paste(plot.dir,'likelihoodslice_SOM.pdf',sep=''),width=4,height=6, colormodel='cmyk', pointsize=11)
par(mfrow=c(2,1), mai=c(.7,.3,.1,.25))
# 240 Myr
plot(mm_example240$co2, mm_example240$fit, type='l', lwd=2, xlim=c(-70,5000), ylim=c(0,8e-4),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
text(3000, 6.8e-4, paste('Likelihood surface slice\n at age =',mm_example240$age,'Myr'))
axis(1, at=seq(0,5000,200), labels=rep('',length(seq(0,5000,200))), col='gray')
axis(1, at=seq(0,5000,1000), labels=c('0','1000','2000','3000','4000','5000'))
arrows(0, 0, 0, 7.7e-4, length=0.08, angle=30, code=2)
mtext(expression(CO[2]~(ppmv)), side=1, line=2.4)
mtext('Density', side=2, line=0.3)
mtext("a.", side=3, line=-0.7, cex=1, adj=-0.06, font=2)
# 140 Myr
plot(mm_example140$co2, mm_example140$fit, type='l', lwd=2, xlim=c(-70,5000), ylim=c(0,8e-4),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
text(3000, 6.8e-4, paste('Likelihood surface slice\n at age =',mm_example140$age,'Myr'))
axis(1, at=seq(0,5000,200), labels=rep('',length(seq(0,5000,200))), col='gray')
axis(1, at=seq(0,5000,1000), labels=c('0','1000','2000','3000','4000','5000'))
arrows(0, 0, 0, 7.7e-4, length=0.08, angle=30, code=2)
mtext(expression(CO[2]~(ppmv)), side=1, line=2.4)
mtext('Density', side=2, line=0.3)
mtext("b.", side=3, line=-0.7, cex=1, adj=-0.06, font=2)
dev.off()



# plot the distribution of modeled CO2 at 240 Myr under likelihood function,
# priors, precalibration and MCMC
plot(mm_modeled$co2, mm_modeled$fit, type='l', lwd=2, lty=3, xlim=c(0,5000),
     xlab='CO2 (ppmv)', ylab='Probability density')
lines(mm_example240$co2, mm_example240$fit, lwd=2)
lines(mm_precal$co2, mm_precal$fit, lwd=2, lty=4)
lines(mm_priors$co2, mm_priors$fit, lwd=2, lty=2)
lines(mm_essgym$co2, mm_essgym$fit, lwd=2, lty=2, col='purple')
legend(2500, 0.0019, c('Likelihood','Priors','Priors, ESS-GYM only','Precal','Posterior'),
       lty=c(1,2,2,4,3), col=c('black','black','purple','black','black'), lwd=2)
##==============================================================================



##==============================================================================
## Other plots showing the supplementary experiments
## with many model configurations
##===================================

TODO

## SOM figure with our ensemble and that of the supplmeental PR2011 experiment
## (only ACT, FERT, GYM, LIFE, deltaT2X and GLAC  parameters calibrated)


## TODO -- need to modify this and the analysis.R script to process the results
## of the many supplemental experiments and get a list of the figures to include
## in SOM


## Log scale (model and points, with CO2 and age error bars, no likelihood surface)
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_PR2011.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', col='white')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.7,.7,.7,.6), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_pr2011[,'q025'],rev(model_quantiles_pr2011[,'q975']))), col=rgb(.94, .6,.6,.6), border=NA)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-457, log10(35), c('95% range, all 69 parameters','95% range, only 6 PR2011 parameters'), pch=c(15,15), col=c('gray','coral'), cex=.75, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


## Log scale (model and points, with CO2 and age error bars, no likelihood surface)
## for experiment with the unimodal likelihood surface, fitted as in
## Park and Royer 2011
pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_PR2011unimodal.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(9000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', col='white')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.7,.7,.7,.7), border=NA)
polygon(-c(time,rev(time)), log10(c(model_quantiles_pr2011uni[,'q025'],rev(model_quantiles_pr2011uni[,'q975']))), col=rgb(.94, .6,.6,.7), border=NA)
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
#for (ii in 1:nrow(data_calib)) {
#    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
#    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)
#}
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-447, log10(32), c('95% range, all 69 parameters','95% range, PR2011 parameter',' and likelihood function'), pch=c(15,15,NA), col=c('gray','coral',NA), cex=.65, bg="white", box.col="white", box.lwd=0)
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()



##==============================================================================
## End
##==============================================================================
