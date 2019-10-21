##==============================================================================
## process_results.R
##
## Read results, compute convergence diagnostics, thin chains. Yields an output
## file to be read into the analysis.R routine.
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

setwd('~/codes/GEOCARB/R')
library('coda')
library('Hmisc')
source("compute_maxlag.R")
source("compute_grdiag.R")
source("pick_chains.R")
chains <- NULL



## [1] #########################################################################
##        dPpPUsOlUsn
appen <- "dPpPUsOlUsn"
datestamp <- "04Oct2019"
source("process_supp_chains.R");


## [2] #########################################################################
##        dPpPUsLlUsn
appen <- "dPpPUsLlUsn"
datestamp <- "04Oct2019"
source("process_supp_chains.R");


## [3] #########################################################################
##        dPpPUsRlUsn
appen <- "dPpPUsRlUsn"
datestamp <- "04Oct2019"
source("process_supp_chains.R");


## [4] #########################################################################
##        dPpPUsRlMsn
appen <- "dPpPUsRlMsn"
datestamp <- "03Oct2019"
source("process_supp_chains.R");


## [5] #########################################################################
##        dFpPUsRlUsn
appen <- "dFpPUsRlUsn"
datestamp <- "03Oct2019"
source("process_supp_chains.R");


## [6] #########################################################################
##        dFpPUsRlMsn
appen <- "dFpPUsRlMsn"
datestamp <- "03Oct2019"
source("process_supp_chains.R");


## [7] #########################################################################
##        dPpAUsRlMsn
## Paper results: use chains 1, 2, 4, 5 (confirmed)
appen <- "dPpAUsRlMsn"
datestamp <- "03Oct2019"
{source("process_supp_chains.R")};


## [8] #########################################################################
##        dFpAUsRlUsn
## Paper results: use chains 1, 3, 4, 5 (confirmed)
appen <- "dFpAUsRlUsn"
datestamp <- "03Oct2019"
{source("process_supp_chains.R")};


## [9] #########################################################################
##        dFpAUsRlMsn
## Paper results: use all chains (confirmed)
appen <- "dFpAUsRlMsn"
datestamp <- "04Oct2019"
{source("process_supp_chains.R")};


## [10] #########################################################################
##        dFpAUsRlMnm
## Paper results: use all chains (confirmed)
appen <- "dFpAUsRlMnm"
datestamp <- "03OCt2019"
{source("process_supp_chains.R")};


################################################################################

# saving for now to work on something else
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_chains <- paste('../output/chains_analysis_',today,'.rds', sep="")
saveRDS(chains, file=filename_chains)

##==============================================================================



##==============================================================================
## END
##==============================================================================
