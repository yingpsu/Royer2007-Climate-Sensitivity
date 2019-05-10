#===============================================================================
# make_initialization_files.R
#
# Create transition covariance matrix and parameters vector for initializing
# Markov chains better.
#
# Tony Wong (anthony.e.wong@colorado.edu)
#===============================================================================

setwd('./codes/GEOCARB/R')

load('../output/geocarb_mcmcoutput_unc_09May2019sn.RData') # based on standard parameters
cov.jump <- amcmc_out1$cov.jump
p0 <- amcmc_out1$samples[nrow(amcmc_out1$samples),]

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
