##==============================================================================
## plotting.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')
plot.dir <- '../figures/'
load('../output/analysis_06Oct2019.RData')

library(Hmisc)

##==============================================================================
# Figure 5 (Methods). Observations and fitted likelihood surface.

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

# only the 5-95% range
# MS FIGURE
pdf(paste(plot.dir,'data_likelihood_withModel_5-95.pdf',sep=''),width=4,height=3,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(likelihood_quantiles[,'50']), type='l', xlim=c(-430,0), ylim=c(0.7,log10(7200)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
grid()
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'05'],rev(likelihood_quantiles[idx_likelihood_ages,'95']))), col='gray70', border=NA)
lines(-time[idx_likelihood_ages], log10(likelihood_quantiles[idx_likelihood_ages,'50']), lwd=2, lty=1)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2, lty=5, col="salmon3")
points(-data_calib$age, log10(data_calib$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=1)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
legend(-428, log10(60), c('Data','Likelihood median','95% range','Model (MLE)'), pch=c(4,NA,15,NA), lty=c(NA,1,NA,2), lwd=c(NA,1.5,NA,1.5), col=c('black','black','gray85','salmon3'), cex=.8, bty='n')
legend(-428, log10(60), c('Data','Likelihood median','95% range','Model (MLE)'), pch=c(NA,NA,NA,NA), lty=c(NA,1,NA,5), lwd=c(NA,1.5,NA,1.5), col=c('black','black','gray85','salmon3'), cex=.8, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()


##==============================================================================
# Figure 2. Posterior model ensemble (gray shaded region denotes 5-95% credible
# range), maximum posterior score simulation (solid bold line) and uncalibrated
# model simulation (dashed line), with proxy data points superimposed (+ markers).

# need a fix for the NAN that results when co2_low for data points is 0
idx_low <- which(data_calib$co2_low == 0)
data_calib$co2_low[idx_low] <- 1

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


## Log scale (model and points, with CO2 and age error bars, no likelihood surface, but with Royer et al 2014 too)
## MS FIGURE 2

pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_2errbars+Royer14.pdf',sep=''),width=4,height=3,colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles_royer[,'q025'],rev(model_quantiles_royer[,'q975']))), col=rgb(.6,.2,.6,.5), border=NA)
lines(-time, log10(model_quantiles_royer[,'co2']), lwd=2, lty=5, col=rgb(.6,.2,.6))
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.2,.6,.6,.5), border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2, col=rgb(.2,.6,.6))
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.25)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.25)
}
#points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='orange', pch=16, cex=1)
#points(-data_calib$age[idx_gastaldo], log10(data_calib$co2[idx_gastaldo]), col='black', pch=16, cex=0.5)
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), las=1)
legend(-452, log10(55), c('Data'), pch=c(16), col=c('black'), pt.cex=0.8, cex=.6, bty='n')
legend(-452, log10(40), c('95% range,\n this work','95% range,\n Royer et al [2014]'), pch=c(15,15), col=c(rgb(.2,.6,.6,.5),rgb(.6,.2,.6,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-452, log10(40), c('95% range,\n this work','95% range,\n Royer et al [2014]'), pch=c('-','-'), col=c(rgb(.2,.6,.6),rgb(.6,.2,.6)), cex=.6, bty='n', y.intersp=1.6)
minor.tick(nx=5, ny=0, tick.ratio=0.5)
dev.off()

## Log scale (model and points, with CO2 and age error bars, no likelihood surface, but with Royer et al 2014 too)
## MS FIGURE SOM

pdf(paste(plot.dir,'model_ensemble_vs_obspts_logscale_2errbars+likelihood.pdf',sep=''),width=4,height=3,colormodel='cmyk', pointsize=11)
par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'maxpost']), type='l', xlim=c(-450,0), ylim=c(0.7,log10(6500)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(likelihood_ages,rev(likelihood_ages)), log10(c(likelihood_quantiles[idx_likelihood_ages,'025'],rev(likelihood_quantiles[idx_likelihood_ages,'975']))), col=rgb(.6,.2,.6,.5), border=NA)
lines(-time[idx_likelihood_ages], log10(likelihood_quantiles[idx_likelihood_ages,'50']), lwd=2, lty=5, col=rgb(.6,.2,.6))
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col=rgb(.2,.6,.6,.5), border=NA)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2, col=rgb(.2,.6,.6))
points(-data_calib$age, log10(data_calib$co2), pch=16, cex=0.4, lwd=.4)
for (ii in 1:nrow(data_calib)) {
    arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
    arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)
}
mtext('Time [Myr ago]', side=1, line=2.1)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)))
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), las=1)
legend(-452, log10(55), c('Data'), pch=c(16), col=c('black'), pt.cex=0.8, cex=.6, bty='n')
legend(-452, log10(40), c('95% range,\n this work','95% range,\n likelihood fcn'), pch=c(15,15), col=c(rgb(.2,.6,.6,.5),rgb(.6,.2,.6,.5)), pt.cex=1.2, cex=.6, bty='n', y.intersp=1.6)
legend(-452, log10(40), c('95% range,\n this work','95% range,\n likelihood fcn'), pch=c('-','-'), col=c(rgb(.2,.6,.6),rgb(.6,.2,.6)), cex=.6, bty='n', y.intersp=1.6)
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
x_pr2011 <- pr2011_icdf(c(.05,.5,.95))
iglac <- match('GLAC',parnames_calib)
x_thisstudy <- quantile(parameters[,ics], c(.05,.5,.95))  # 3.148759 4.203719 5.621429
x_glac <- quantile(parameters[,ics]*parameters[,iglac], c(.05,.5,.95))  # 5.678216  9.289884 14.268706
x_norm <- quantile(parameters_nm[,ics], c(.05,.5,.95))  # 3.273916 4.396285 5.846267
x_ktc2017 <- c(3.7, 5.6, 7.5)

offset <- 0.1

## MS FIGURE

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
y0 <- 0.7*offset; arrows(x_thisstudy[1], y0, x_thisstudy[3], y0, lwd=1.5, length=0.04, angle=90, code=3, col="steelblue"); points(x_thisstudy[2], y0, pch=16, col="steelblue")
#y1 <- 0.35*offset; arrows(x_royer2007[1], y1, x_royer2007[3], y1, lwd=1.5, length=0.04, angle=90, code=3); points(x_royer2007[2], y1, pch=15)
y1 <- 0.3*offset; arrows(x_pr2011[1], y1, x_pr2011[3], y1, lwd=1.5, length=0.04, angle=90, code=3); points(x_pr2011[2], y1, pch=15)
#y2 <- 0.08; arrows(x_ktc2017[1], y2, x_ktc2017[3], y2, length=0.04, angle=90, code=3); points(x_ktc2017[2], y2, pch=17)
legend(5.1,1.02, c('5-95% range, PR2011','PR2011','5-95% range, this study','Posterior, this study','Prior, both studies'),
       pch=c(15,NA,16,NA,NA), lty=c(1,2,1,1,3), col=c("black","black","steelblue","steelblue","steelblue"), bty='n', lwd=1.7, cex=0.9)
dev.off()


##==============================================================================
# MS FIGURE SOM. Posterior probability density for Earth system sensitivity
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
print(paste('PR2011 dT2X 5-50-95% quantiles =',pr2011_icdf(c(.05,.5,.95))))
print(paste('PR2011 dT2X 16-50-84% quantiles =',pr2011_icdf(c(.16,.5,.84))))

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

## MS FIGURE SOM

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
plot(mm_example50$co2, mm_example50$fit, type='l', lwd=2, xlim=c(-70,5000), ylim=c(0,2.5e-3),
     xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n', axes=FALSE)
text(3000, 2e-3, paste('Likelihood surface slice\n at age =',mm_example50$age,'Myr'))
axis(1, at=seq(0,5000,200), labels=rep('',length(seq(0,5000,200))), col='gray')
axis(1, at=seq(0,5000,1000), labels=c('0','1000','2000','3000','4000','5000'))
arrows(0, 0, 0, 2.4e-3, length=0.08, angle=30, code=2)
mtext(expression(CO[2]~(ppmv)), side=1, line=2.4)
mtext('Density', side=2, line=0.3)
mtext("b.", side=3, line=-0.7, cex=1, adj=-0.06, font=2)
dev.off()

# not included, but nice to have:
if(FALSE){
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
}
##==============================================================================



##==============================================================================
## MS FIGURE SOM -- supplementary experiment model vs Royer 2014 comparisons
##===================================

# SOM FIGURE - 5 row x 2 col figure of each experiment relative to the original Royer
# et al 2014 calibration results. Include the data set for each experiment.

pdf('../figures/model_experiments_vs_r2014.pdf',width=7,height=9,colormodel='cmyk', pointsize=9)
par(mfrow=c(5,2))
par(cex = 0.85)
par(mar = c(0, 0, 2.5, 0), oma = c(3.5, 5, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
panel <- 1
panel_label <- c("a.","b.","c.","d.","e.","f.","g.","h.","i.","j.")
for (ee in names(chains)) {
    chainname <- ee
    if (chainname == "x_pr2011") {table_experiments[row,] <- c("PR2011", "PR2011", "PR2011", "unimodal", "log-normal")
    } else {
        if (substr(chainname, 2, 2)=="P") {data_choice <- "PR2011"} else {data_choice <- "F2017"}
        if (substr(chainname, 4, 5)=="PU") {param_choice <- "PR2011"} else {param_choice <- "all"}
        if (substr(chainname, 7, 7)=="R") {fSR_choice <- "DT2019"} else if (substr(chainname, 7, 7)=="L") {fSR_choice <- "L2018"} else {fSR_choice <- "PR2011"}
        if (substr(chainname, 9, 9)=="U") {likelihood_choice <- "unimodal"} else {likelihood_choice <- "mixture"}
        if (substr(chainname, 10, 11)=="sn") {kernel_choice <- "skew-normal"} else {kernel_choice <- "normal"}
    }
    title_tag <- paste(data_choice,"data,",param_choice,"parameters,",fSR_choice,"fSR,",likelihood_choice,"likelihood,",kernel_choice,"uncertainties")
    title_tag1 <- paste(data_choice,"data,",param_choice,"parameters,",fSR_choice,"fSR")
    title_tag2 <- paste(likelihood_choice,"likelihood,",kernel_choice,"uncertainties")

    if (substr(ee, 2,2)=="P") {data_calib <- data_calib_pr2011
    } else {data_calib <- data_calib_f2017}
    plot(-time, log10(model_experiment_quantiles[[ee]][,'0.5']), type='l', xlim=c(-450,0), ylim=c(2,log10(6000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
    points(-data_calib$age, log10(data_calib$co2), pch=4, cex=0.4, lwd=.4)
    for (ii in 1:nrow(data_calib)) {
        arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
        if(substr(ee, 2,2)=="F") {arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)}
    }
    polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles$r2014[,'0.025'],rev(model_experiment_quantiles$r2014[,'0.975']))), col=rgb(.7,.2,.4,.5), border=NA)
    lines(-time, log10(model_experiment_quantiles$r2014[,'0.5']), lwd=1.5, lty=5, col=rgb(.6,.2,.6))
    polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles[[ee]][,'0.025'],rev(model_experiment_quantiles[[ee]][,'0.975']))), col=rgb(.2,.6,.6,.5), border=NA)
    lines(-time, log10(model_experiment_quantiles[[ee]][,'0.5']), lwd=1.5, col=rgb(.2,.6,.75))
    #mtext(ee, side=3)
    mtext(panel_label[panel], side=3, adj=0.01, font=2)
    mtext(title_tag1, side=3, line=0.8)
    mtext(title_tag2, side=3, line=0)
    if (panel %in% c(9,10)) {
        axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'))
        mtext('Time [Myr ago]', side=1, line=2.1)
    }
    ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
    if (panel %in% c(1,3,5,7,9)) {
        axis(2, at=ticks, labels=rep('',length(ticks)))
        axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), las=1)
        mtext(expression('CO'[2]*' [ppmv]'), side=2, line=3)
    }
    if (panel==1) {
        legend(-390, log10(6600), c('95% range (Royer [2014])', '95% range (this work)'), pch=c(15,15), col=c(rgb(.7,.2,.4,.5), rgb(.2,.6,.6,.5)), cex=.8, bty='n')
        legend(-180, log10(6600), c('',''), pch=c(NA,4), col=c('black','black'), cex=.8, bty='n')
        legend(-180, log10(6600), c('Median','Data'), pch=c('-',NA), col=c('black','black'), cex=.8, bty='n')
    }
    minor.tick(nx=5, ny=0, tick.ratio=0.5)
    panel <- panel+1
}
dev.off()

##==============================================================================



##==============================================================================
## MS FIGURE SOM -- quantiles of deltaT2X as sample size increases
##===========

col_x <- c(8, 10.4, 13.3, 15.8, 18.5)
col_names <- c("Data", "Parameters", "fSR", "Likelihood", "Uncertainties")

table_experiments <- mat.or.vec(nr=nrow(ess$ess), length(col_names))
colnames(table_experiments) <- col_names
for (row in 1:nrow(ess$ess)) {
    chainname <- rownames(ess$ess)[row]
    if (chainname == "x_pr2011") {table_experiments[row,] <- c("PR2011", "PR2011", "PR2011", "unimodal", "log-normal")
    } else {
        if (substr(chainname, 2, 2)=="P") {table_experiments[row, "Data"] <- "PR2011"} else {table_experiments[row, "Data"] <- "F2017"}
        if (substr(chainname, 4, 5)=="PU") {table_experiments[row, "Parameters"] <- "PR2011"} else {table_experiments[row, "Parameters"] <- "all"}
        if (substr(chainname, 7, 7)=="R") {table_experiments[row, "fSR"] <- "DT2019"} else if (substr(chainname, 7, 7)=="L") {table_experiments[row, "fSR"] <- "L2018"} else {table_experiments[row, "fSR"] <- "PR2011"}
        if (substr(chainname, 9, 9)=="U") {table_experiments[row, "Likelihood"] <- "unimodal"} else {table_experiments[row, "Likelihood"] <- "mixture"}
        if (substr(chainname, 10, 11)=="sn") {table_experiments[row, "Uncertainties"] <- "skew-normal"} else {table_experiments[row, "Uncertainties"] <- "normal"}
    }
}

pdf(file='../figures/boxplot_ess.pdf', width=8, height=3.7, colormodel="cmyk", pointsize=11)
offset <- 0.045
yhgt <- offset*.75
experiment <- rownames(ess$ess)[1]
par(mfrow=c(1,1), mai=c(.7,.2,.2,.2))
plot(ess$ess[experiment,"0.5"], yhgt, xlim=c(0,22), ylim=c(0, .58), pch=16,
     xaxs='i', yaxs='i', yaxt='n', ylab='', xlab='', axes=FALSE)
grid()
axis(1, at=seq(0,9,2)); axis(1, at=seq(0,8,1), labels=rep("", 9), lwd=0.25)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.2, adj=0.18)
for (experiment in rownames(ess$ess)) {
    row <- which(rownames(ess$ess)==experiment)
    points(ess$ess[experiment,"0.5"], yhgt, pch=16)
    arrows(x0=ess$ess[experiment,"0.05"], x1=ess$ess[experiment,"0.95"], y0=yhgt, y1=yhgt, angle=90, length=0.05, code=3)
    #text(8, yhgt, experiment, pos=4)
    #text(8, yhgt, table_experiments[row,], pos=4)
    for (i in 1:length(col_x)) {text(col_x[i], yhgt, table_experiments[row, i], pos=4)}
    yhgt <- yhgt + offset
}
yhgt <- yhgt + offset*0.25
for (i in 1:length(col_x)) {text(col_x[i], yhgt, colnames(table_experiments)[i], pos=4)}
lines(c(8, 22), c(yhgt-0.4*offset, yhgt-0.4*offset))
dev.off()

##==============================================================================



##==============================================================================
## MS FIGURE SOM -- quantiles of deltaT2X as sample size increases
##===========

pdf(file='../figures/deltaT2X_quantiles.pdf', width=3.5, height=3, colormodel="cmyk", pointsize=11)
par(mfrow=c(1,1), mai=c(.7,.7,.2,.3))
plot(sample_sizes, sample_quantiles[,1], type='l', ylim=c(0,8), lty=2, xlab='', ylab='')
lines(sample_sizes, sample_quantiles[,3], lty=2)
lines(sample_sizes, sample_quantiles[,2])
grid()
mtext("Sample size", side=1, line=2.2)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=2, line=2.2)
dev.off()

##==============================================================================


##==============================================================================
## End
##==============================================================================
