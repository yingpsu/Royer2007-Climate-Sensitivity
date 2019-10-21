##==============================================================================
## likelihood_surface_quantiles.R
##
## Assuming fit_likelihood_surface has been run previously, compute the
## quantiles for each time slice from the mixture model likelihood (or whichever
## likelihood was fit).
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

##==============================================================================
## Sample and fit KDEs to each time step
##======================================

library(Bolstad)

likelihood_quantiles <- array(NA, c(n_time, 10))
colnames(likelihood_quantiles) <- c('025','05','17','25','50','75','83','95','975','Max')
quantiles_to_compute <- c(.025,.05,.17,.25,.5,.75,.83,.95,.975)
x_co2 <- seq(from=0, to=10000, by=1)

for (tt in 1:n_time) {
  if (!is.null(likelihood_fit[[tt]])) {
    tmp <- sintegral(x_co2, likelihood_fit[[tt]])
    for (qq in 1:length(quantiles_to_compute)) {
      likelihood_quantiles[tt,qq] <- tmp$cdf$x[which.min(abs(tmp$cdf$y-quantiles_to_compute[qq]))]
    }
    # max likelihood
    lhood <- likelihood_fit[[tt]](x_co2)
    likelihood_quantiles[tt,10] <- x_co2[which.max(lhood)]
  }
}

idx_likelihood_ages <- which(!is.na(likelihood_quantiles[,1]))
likelihood_ages <- time[idx_likelihood_ages]

##==============================================================================


##==============================================================================
## End
##==============================================================================
