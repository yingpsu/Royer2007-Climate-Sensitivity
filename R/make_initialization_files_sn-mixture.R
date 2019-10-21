##==============================================================================
## make_initialization_files_sn-mixture.R
##
## Create transition covariance matrix and parameters vector for initializing
## Markov chains better.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

setwd('~/codes/GEOCARB/R')

appen <- 'sn-mix'
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
