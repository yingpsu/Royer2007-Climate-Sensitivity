#===============================================================================
# preliminary_scratch.R
#===============================================================================

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
plot.dir <- '~/Google\ Drive/SCRiM/projects/geocarb/figures/'

#===============================================================================
# results calibrating all of the 56 constant parameters
#

# read some prelminary results in
load('../output/GEOCARB_MCMC_allConst-alkenones_g60_09Nov2017.RData');  chain1_alk = chain1
load('../output/GEOCARB_MCMC_allConst-paleosols_g60_09Nov2017.RData');  chain1_pal = chain1
load('../output/GEOCARB_MCMC_allConst-liverworts_g60_09Nov2017.RData'); chain1_liv = chain1
load('../output/GEOCARB_MCMC_allConst-boron_g60_09Nov2017.RData');      chain1_bor = chain1
load('../output/GEOCARB_MCMC_allConst-stomata_g60_09Nov2017.RData');    chain1_sto = chain1
load('../output/GEOCARB_MCMC_allConst-allData_g60_09Nov2017.RData');    chain1_all = chain1
load('../output/GEOCARB_MCMC_allConst-paleosols_g50_09Nov2017.RData');  chain1_pal_g50 = chain1

# which index is deltaT2X
ics <- which(parnames_calib=='deltaT2X')

pdf(paste(plot.dir,'historyPlot_allConst.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.4,.6,.6,.2))
  plot(chain1_pal[,ics], type='l', xlab='', ylab='CS')
  mtext('Paleosols', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
  mtext('Calibrating all constant parameters', side=3, line=2.8, adj=0, cex=.9)
par(mai=c(.4,.6,.6,.2))
  plot(chain1_sto[,ics], type='l', xlab='', ylab='CS')
  mtext('Stomata', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_bor[,ics], type='l', xlab='', ylab='CS')
  mtext('Boron', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_alk[,ics], type='l', xlab='', ylab='CS')
  mtext('Alkenones', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_liv[,ics], type='l', xlab='', ylab='CS')
  mtext('Liverworts', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_all[,ics], type='l', xlab='', ylab='CS')
  mtext('All data', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
dev.off()

# plot only the paleosols using slower/faster rates of adaptation
pdf(paste(plot.dir,'historyPlot_allConst_adaptationRate.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.4,.6,.6,.2))
  plot(chain1_pal[,ics], type='l', xlab='', ylab='CS')
  mtext('Paleosols, slower adaptation', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
  mtext('Calibrating all constant parameters', side=3, line=2.8, adj=0, cex=.9)
par(mai=c(.4,.6,.6,.2))
  plot(chain1_pal_g50[,ics], type='l', xlab='', ylab='CS')
  mtext('Paleosols, faster adaptation', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
dev.off()
#===============================================================================


#===============================================================================
# results calibrating just a handful of constant parameters
#

# read some prelminary results in
load('../output/GEOCARB_MCMC_fewConst-alkenones_g50_10Nov2017.RData');  chain1_alk = chain1
load('../output/GEOCARB_MCMC_fewConst-paleosols_g50_10Nov2017.RData');  chain1_pal = chain1
load('../output/GEOCARB_MCMC_fewConst-liverworts_g50_10Nov2017.RData'); chain1_liv = chain1
load('../output/GEOCARB_MCMC_fewConst-boron_g50_10Nov2017.RData');      chain1_bor = chain1
load('../output/GEOCARB_MCMC_fewConst-stomata_g50_10Nov2017.RData');    chain1_sto = chain1
load('../output/GEOCARB_MCMC_fewConst-allData_g50_10Nov2017.RData');    chain1_all = chain1

ics <- which(parnames_calib=='deltaT2X')

# plot em
pdf(paste(plot.dir,'historyPlot_fewConst.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.4,.6,.6,.2))
  plot(chain1_pal[,ics], type='l', xlab='', ylab='CS')
  mtext('Paleosols', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
  mtext('Calibrating a few constant parameters', side=3, line=2.8, adj=0, cex=.9)
par(mai=c(.4,.6,.6,.2))
  plot(chain1_sto[,ics], type='l', xlab='', ylab='CS')
  mtext('Stomata', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_bor[,ics], type='l', xlab='', ylab='CS')
  mtext('Boron', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_alk[,ics], type='l', xlab='', ylab='CS')
  mtext('Alkenones', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_liv[,ics], type='l', xlab='', ylab='CS')
  mtext('Liverworts', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_all[,ics], type='l', xlab='', ylab='CS')
  mtext('All data', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
dev.off()
#===============================================================================


#===============================================================================
# results calibrating just a handful of constant parameters AND a handful of
# time-varying ones
#

# read some prelminary results in
load('../output/GEOCARB_MCMC_fewBoth-alkenones_g50_11Nov2017.RData');  chain1_alk = chain1
load('../output/GEOCARB_MCMC_fewBoth-paleosols_g50_11Nov2017.RData');  chain1_pal = chain1
load('../output/GEOCARB_MCMC_fewBoth-liverworts_g50_11Nov2017.RData'); chain1_liv = chain1
load('../output/GEOCARB_MCMC_fewBoth-boron_g50_11Nov2017.RData');      chain1_bor = chain1
load('../output/GEOCARB_MCMC_fewBoth-stomata_g50_11Nov2017.RData');    chain1_sto = chain1
load('../output/GEOCARB_MCMC_fewBoth-allData_g50_11Nov2017.RData');    chain1_all = chain1

ics <- which(parnames_calib=='deltaT2X')

# plot em
pdf(paste(plot.dir,'historyPlot_fewBoth.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(3,2), mai=c(.4,.6,.6,.2))
  plot(chain1_pal[,ics], type='l', xlab='', ylab='CS')
  mtext('Paleosols', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
  mtext('Calibrating a few constant and a few time-varying parameters', side=3, line=2.8, adj=0, cex=.9)
par(mai=c(.4,.6,.6,.2))
  plot(chain1_sto[,ics], type='l', xlab='', ylab='CS')
  mtext('Stomata', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_bor[,ics], type='l', xlab='', ylab='CS')
  mtext('Boron', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.45,.6,.55,.2))
  plot(chain1_alk[,ics], type='l', xlab='', ylab='CS')
  mtext('Alkenones', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_liv[,ics], type='l', xlab='', ylab='CS')
  mtext('Liverworts', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
par(mai=c(.5,.6,.5,.2))
  plot(chain1_all[,ics], type='l', xlab='', ylab='CS')
  mtext('All data', side=3, line=.8, cex=.8); mtext('Iteration', side=1, line=2.2, cex=0.75)
dev.off()
#===============================================================================
