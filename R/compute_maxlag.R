##==============================================================================
## compute_maxlag.R
##
## Computes the lag at which the Markov chains stored in the input `chains` are
## relatively non-autocorrelated. `chains` is assumed to be a list object with
## each list element corresponding to a matrix of Markov chain iterates, where
## each row in the matrix is a concommitant parameter sample, and the number of
## rows = the number of Markov chain iterations. `chains` can also be an input
## matrix where nrow = the number of iterations and ncol = the number of
## parameters.
## A correlation threshold of 5% is hard-coded (`cmax`, below). The maximum lag
## considered for each parameter for each chain starts out at `lmax0`, and is
## increased (by 200, by default) if the parameter chain does not have any lag
## below `lmax` at which the autocorrelation is below `cmax`.
## This is then reset to `lmax0` for the next set of Markov chains (in the input
## list `chains`). You do NOT want to continue at the previous chain's final
## value for `lmax` because you will end up computing a LOT of extra
## autocorrelation values, and the function will take notably more time to run.
## The maximum lag from each chain's max lag may be used to thin the chains to
## obtain (relatively) independent samples from their stationary distribution.
##
## Tony Wong (aewsma@rit.edu)
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

compute_maxlag <- function(chains) {
  lmax0 <- 1000
  cmax <- 0.05

  if (is.list(chains)) {
    maxlags <- rep(0, length(chains))
    for (m in 1:length(chains)) {
        lmax <- lmax0
        for (p in 1:ncol(chains[[m]])) {
            acf_tmp <- acf(chains[[m]][,p], lag.max=lmax, plot=FALSE)
            idx_low <- which(acf_tmp$acf < cmax)
            while (length(idx_low)==0) {
              lmax <- lmax + 200
              acf_tmp <- acf(chains[[m]][,p], lag.max=lmax, plot=FALSE)
              idx_low <- which(acf_tmp$acf < cmax)
            }
            new <- acf_tmp$lag[idx_low[1]]
            if (maxlags[m] < new) {
                print(paste(m,p,"Updating maxlags[m] to",new))
                maxlags[m] <- new
            }
        }
    }
    return(maxlags)
  } else {
    maxlag <- 0
    lmax <- lmax0
    for (p in 1:ncol(chains)) {
        acf_tmp <- acf(chains[,p], lag.max=lmax, plot=FALSE)
        idx_low <- which(acf_tmp$acf < cmax)
        while (length(idx_low)==0) {
          lmax <- lmax + 200
          acf_tmp <- acf(chains[,p], lag.max=lmax, plot=FALSE)
          idx_low <- which(acf_tmp$acf < cmax)
        }
        new <- acf_tmp$lag[idx_low[1]]
        if (maxlag < new) {
            print(paste(p,"Updating maxlag to",new))
            maxlag <- new
        }
    }
    return(maxlag)
  }
}

##==============================================================================
## End
##==============================================================================
