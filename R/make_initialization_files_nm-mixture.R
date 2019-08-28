##==============================================================================
## make_initialization_files_nm-mixture.R
##
## Create transition covariance matrix and parameters vector for initializing
## Markov chains better. This version is for the supplemental experiment using
## a Gaussian mixture model (normal).
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

setwd('~/codes/GEOCARB/R')

appen <- 'nm-mix'
load('../output/geocarb_mcmcoutput_mix_06Jul2019sn-mix.RData') # based on standard parameters, an initial 1e6 run
parameters <- amcmc_out1$samples
load('../output/geocarb_mcmcoutput_mix_07Jul2019sn-mix.RData') # based on previous run, another 2e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_mix_08Jul2019sn-mix.RData') # based on previous run, another 2e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_mix_12Jul2019sn-mix.RData') # based on previous run, another 2e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_mix_15Jul2019sn-mix.RData') # based on previous run, another 2e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_mix_20Jul2019sn-mix.RData') # based on previous run, another 2e6
parameters <- rbind(parameters, amcmc_out1$samples)
load('../output/geocarb_mcmcoutput_mix_10Aug2019sn-mix.RData') # based on previous run, another 1e6
parameters <- rbind(parameters, amcmc_out1$samples)

# create initialization files from the end of the last chain read

cov.jump <- amcmc_out1$cov.jump
p0 <- amcmc_out1$samples[nrow(amcmc_out1$samples),]

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
