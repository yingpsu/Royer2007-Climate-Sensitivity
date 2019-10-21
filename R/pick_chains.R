##==============================================================================
## pick_chains.R
##
## To be run after plotting all of the parallel Markov chains for an experiment.
## Gives the user the option of which chains to keep for analysis.
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

pick_chains <- function(n) {
    use_chains <- rep(TRUE, n)
    {for (m in 1:n) {use_chains[m] <- as.logical(readline(prompt=paste("Use chain",m,"? Enter TRUE or FALSE. ")))}};
    return(use_chains)
}

##==============================================================================
## End
##==============================================================================
