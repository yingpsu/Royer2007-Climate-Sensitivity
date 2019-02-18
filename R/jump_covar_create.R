#===============================================================================
# jump_covar_create.R
#
# Create transition covariance matrix for initializing Markov chains better.
#
# Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

setwd('./codes/GEOCARB/R')

load('../output/GEOCARB_MCMC_unc-sd10_29Jan2019sn.RData')

appen <- 'unc-sd10'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
output_dir <- '../output/'

# write transition covariance matrix initialization file
filename_cov <- paste(output_dir,'covar_init_',appen,'_',today,'.rds', sep='')
saveRDS(object=amcmc_out1$cov.jump, file=filename_cov)

# write parameters initialization file
filename_par <- paste(output_dir,'param_init_',appen,'_',today,'.rds', sep='')
saveRDS(object=amcmc_out1$samples[nrow(amcmc_out1$samples),], file=filename_par)

#===============================================================================
# End
#===============================================================================
