#===============================================================================
# jump_covar_create.R
#
# Create transition covariance matrix for initializing Markov chains better.
#
# Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

setwd('./codes/GEOCARB/R')

#load('../output/GEOCARB_MCMC_unc-sd10_29Jan2019sn.RData')  # first time
load('../output/geocarb_mcmcoutput_unc_04Apr2019sn.RData') # third time, based on second
cov.jump <- amcmc_out1$cov.jump
p0 <- amcmc_out1$samples[nrow(amcmc_out1$samples),]

#load('../output/geocarb_mcmcoutput_unc_25Mar2019sn.RData') # second time, based on first
#cov.jump <- amcmc_par1[[1]]cov.jump
#p0 <- amcmc_par1[[1]]$samples[nrow(amcmc_par1[[1]]$samples),]

appen <- 'unc-sd10'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
output_dir <- '../output/'

# write transition covariance matrix initialization file
filename_cov <- paste(output_dir,'covar_init_',appen,'_',today,'.rds', sep='')
saveRDS(object=cov.jump, file=filename_cov)

# write parameters initialization file
filename_par <- paste(output_dir,'param_init_',appen,'_',today,'.rds', sep='')
saveRDS(object=p0, file=filename_par)

#===============================================================================
# End
#===============================================================================
