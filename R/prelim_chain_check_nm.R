##==============================================================================
## prelim_chain_check_nm.R
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

setwd('~/codes/GEOCARB/R')

## load the MCMC results and paste together (they are continutations of one another)
load('../output/geocarb_mcmcoutput_unc_28May2019sn.RData') # based on standard parameters, from 3e6 run
parameters <- amcmc_out1$samples
load('../output/geocarb_mcmcoutput_unc_30May2019sn.RData') # based on previous run, another 3e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_unc_02Jun2019sn.RData') # based on previous run, another 3e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_unc_04Jun2019sn.RData') # based on previous run, another 3e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_unc_12Jun2019sn.RData') # based on previous run, another 3e6
parameters <- rbind(parameters, amcmc_out1$samples)
# here is where we shift to a normal distribution
load('../output/geocarb_mcmcoutput_unc_20Jun2019nm.RData') # based on previous run, another 3e6
parameters <- rbind(parameters, amcmc_out1$samples)

# get and set up the parameters
USE_LENTON_FSR <- FALSE
USE_DT2019_FSR <- TRUE
filename.calibinput <- "../input_data/GEOCARB_input_summaries_calib_unc.csv"
source('parameterSetup_tvq.R')

##==============================================================================
## Start off by just straight-up chopping off the first 1e7

parameters0 <- parameters
parameters <- parameters0[(1e7+1):nrow(parameters0),]

##==============================================================================
## Burn in

library(coda)
library(Hmisc)
library(sn)

n_parameters <- ncol(parameters)
niter.test <- c(0)
gr.test <- rep(0, length(niter.test))

n_sample <- nrow(parameters)
last_bit <- 3e6
n_half <- 0.5*last_bit

#chain1 <- parameters[1:n_half,]
#chain2 <- parameters[(n_half+1):n_sample,]
chain1 <- parameters[(n_sample-last_bit+1):(n_sample-last_bit+n_half),]
chain2 <- parameters[(n_sample-last_bit+n_half+1):n_sample,]
niter_mcmc <- nrow(chain1)

string.mcmc.list <- 'mcmc1'
for (m in 2:2) {
  string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
}
for (i in 1:length(niter.test)) {
  for (m in 1:2) {
    eval(parse(text=paste('mcmc',m,' <- as.mcmc(chain',m,'[(niter.test[i]+1):niter_mcmc,])', sep='')))
  }
  eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))
  gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
}

# it's a TINY number!  with results through 12Jun2019, and using the last
# 2e6 samples cut in half, get:
#> gr.test
#[1] 1.016068

chains_burned <- vector('list', 2)
chains_burned[[1]] <- chain1
chains_burned[[2]] <- chain2

##==============================================================================
## Thinning
##=========

lmax <- 3000
cmax <- 0.05
maxlag <- 0

for (m in 1:2) {
  for (p in 1:n_parameters) {
    acf_tmp <- acf(chains_burned[[m]][,p], lag.max=lmax, plot=FALSE)
    new <- acf_tmp$lag[which(acf_tmp$acf < cmax)[1]]
    if (maxlag < new) {
      print(paste(m,p,"Updating maxlag to",new))
      maxlag <- new
    }
  }
}

chains_burned_thinned <- chains_burned # initialize
for (m in 1:2) {
  chains_burned_thinned[[m]] <- chains_burned[[m]][seq(from=1, to=nrow(chains_burned[[m]]), by=maxlag),]
}

parameters_posterior <- chains_burned_thinned[[1]]
for (m in 2:2) {
  parameters_posterior <- rbind(parameters_posterior, chains_burned_thinned[[m]])
}

##==============================================================================
## Run the model ensemble

dist <- 'sn'
DO_SAMPLE_TVQ <- TRUE
USE_LENTON_FSR <- FALSE
USE_DT2019_FSR <- TRUE
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_unc.csv'

source('GEOCARB-2014_parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')
#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `var` uncertainty statistical parameter
source('GEOCARB_fit_likelihood_surface.R')
#source('likelihood_surface_quantiles.R')

# Get model parameter prior distribution bounds
names <- as.character(input$parameter)
bound_lower <- rep(NA, length(names))
bound_upper <- rep(NA, length(names))

ind_neg_inf <- which(input[,'lower_limit']=='_inf')
bound_lower[ind_neg_inf] <- -Inf
bound_lower[setdiff(1:length(names), ind_neg_inf)] <- as.numeric(as.character(input$lower_limit[setdiff(1:length(names), ind_neg_inf)]))
bound_upper <- input$upper_limit

bounds <- cbind(bound_lower, bound_upper)
rownames(bounds) <- as.character(input$parameter)

# only actually need the calibration parameters' bounds, so reformat the bounds
# array to match the vector of calibration parameters
bounds_calib <- mat.or.vec(nr=length(parnames_calib), nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:length(parnames_calib)) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

rm(list=c('bound_lower','bound_upper','bounds'))

## Actually run the ensemble
model_out <- sapply(X=1:nrow(parameters_posterior),
                    FUN=function(k){model_forMCMC(par_calib=parameters_posterior[k,],
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
n_time <- nrow(model_out)
n_ensemble <- nrow(parameters_posterior)
n_parameter <- ncol(parameters_posterior)

quantiles_i_want <- c(0,0.005,.025,.05,.5,.95,.975,0.995,1)
model_quantiles <- mat.or.vec(nr=n_time, nc=(length(quantiles_i_want)+1))
colnames(model_quantiles) <- c('q000','q005','q025','q05','q50','q95','q975','q995','q100','maxpost')
for (t in 1:n_time) {
  #model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_out_noisy[t,], quantiles_i_want)
  model_quantiles[t,1:length(quantiles_i_want)] <- quantile(model_out[t,], quantiles_i_want)
}

model_ref <- model_forMCMC(par_calib=par_calib0,
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
                           par_time_stdev=par_time_stdev)[,'co2']

#par(mfrow=c(1,1), mai=c(.65,.9,.15,.15))
plot(-time, log10(model_quantiles[,'q50']), type='l', xlim=c(-450,0), ylim=c(2,log10(6900)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q000'],rev(model_quantiles[,'q100']))), col='gray', border=NA)
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q005'],rev(model_quantiles[,'q995']))), col='gray', border=NA)
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q025'],rev(model_quantiles[,'q975']))), col='gray', border=NA)
#polygon(-c(time,rev(time)), log10(c(model_quantiles[,'q05'],rev(model_quantiles[,'q95']))), col='gray', border=NA)
#lines(-time, log10(likelihood_quantiles[,'50']), lwd=2, lty=2)
lines(-time, log10(model_quantiles[,'maxpost']), lwd=2)
#lines(-time, model_ref, lwd=2, lty=2)
points(-data_calib$age, log10(data_calib$co2), pch='x', cex=0.65)
mtext('Time [Myr ago]', side=1, line=2.1, cex=.9)
mtext(expression('CO'[2]*' concentration [ppmv]'), side=2, line=3.2, cex=.9)
axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1)
ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1)
axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1, las=1)
legend(-400, log10(8000), c('95% credible range','Max posterior','Data'), pch=c(15,NA,4), col=c('gray','black','black'), cex=.8, bty='n')
legend(-400, log10(8000), c('95% credible range','Max posterior','Data'), pch=c(NA,'-',NA), col=c('gray','black','black'), cex=.8, bty='n')
minor.tick(nx=5, ny=0, tick.ratio=0.5)

##==============================================================================
## End
##==============================================================================
