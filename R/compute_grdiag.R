##==============================================================================
## compute_grdiag.R
##
## Compute the Gelman and Rubin (1992) potential scale reduction factor for a
## set of Markov chains.
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

compute_grdiag <- function(chains, ibeg) {

  string.mcmc.list <- 'mcmc1'
  for (m in 2:length(chains)) {
      string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
  }
  for (m in 1:length(chains)) {
      # convert each of the chains into mcmc object
      eval(parse(text=paste('mcmc',m,' <- as.mcmc(chains[[m]][(ibeg+1):nrow(chains[[m]]),])', sep='')))
  }
  eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))
  gr.test <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
  return(gr.test)
}
