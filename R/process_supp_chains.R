##==============================================================================
## process_supp_chains.R
##
## Performs the actual burn-in and thinning calculations for each supplemental
## experiment set of Markov chains.
##
## Called from supplemental_experiments_processing.R
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

load(paste('../output/geocarb_mcmcoutput_',appen,'_',datestamp,'.RData', sep=''))
niter_mcmc <- nrow(amcmc_par1[[1]]$samples)
n_parameters <- ncol(amcmc_par1[[1]]$samples)
n_node000 <- length(amcmc_par1)

# if only 7 parameters, just use all of the chains
# if using all 69 parameters, thinning takes a long time, so be more strategic
if (n_parameters==7) {
  idx_chains <- NULL
} else {
  # make a plot and offer user a chance to say which chains they want
  if (n_parameters==7) {idx_ess <- 5
  } else {idx_ess <- 10}
  par(mfrow=c(n_node000,1)); for (m in 1:n_node000) {plot(amcmc_par1[[m]]$samples[,idx_ess], type='l')}
  {use_chains <- pick_chains(n_node000)};
  idx_chains <- which(use_chains)
}

# burn-in
ibeg <- 0
if (!is.null(idx_chains)) {
  chains_burned <- vector('list',length(idx_chains)); for (m in 1:length(idx_chains)) {chains_burned[[m]] <- amcmc_par1[[idx_chains[m]]]$samples[(ibeg+1):niter_mcmc,]}
} else {
  chains_burned <- vector('list',n_node000); for (m in 1:n_node000) {chains_burned[[m]] <- amcmc_par1[[m]]$samples[(ibeg+1):niter_mcmc,]}
}
gr.test <- compute_grdiag(chains_burned, ibeg)
print(paste(appen,"-- PSRF =",gr.test))
if(gr.test < 1.1) {print(paste(appen,"-- CONVERGED."))} else {print(paste("WARNING:",appen,"-- NOT CONVERGED."))}

# thinning
maxlags <- compute_maxlag(chains_burned)

# process, and report out
chains[[appen]] <- chains_burned[[1]][seq(from=1, to=niter_mcmc, by=maxlags[1]),]
for (m in 2:length(chains_burned)) {chains[[appen]] <- rbind(chains[[appen]], chains_burned[[m]][seq(from=1, to=niter_mcmc, by=maxlags[m]),])}
print(paste(appen,"-- num_samples =",nrow(chains[[appen]])))

##==============================================================================
## End
##==============================================================================
