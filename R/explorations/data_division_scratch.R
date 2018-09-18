#===============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')

.niter_mcmc <- 2e4   # number of MCMC iterations per node (Markov chain length)
.n_node <- 1         # number of CPUs to use
.n_chain <- 1          # number of parallel MCMC chains, per shard (subsample)
#.n_data <- 50       # number of data points to use in each shard
.n_shard <- 30      # number of data subsamples to use and recombine with consensus MC
gamma_mcmc <- 0.66

# if both are FALSE, then shards are randomly sampled
# if both are TRUE, then we break by data type first, then by time, assuming
# that .n_shard is a multiple of 5 (b/c 5 data types)
break_time <- TRUE  # break shards across time
break_type <- TRUE  # break shards across proxy data types

#appen <- 'sig18+GLAC+LIFE'
appen <- 'sig18'
#appen <- 'all'
appen2 <- 'even'
output_dir <- '../output/'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
co2_uncertainty_cutoff <- 20

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

#filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_GAMMA-co2_31Jul2018.csv'
#filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_LN-co2_31Jul2018.csv'
#filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_06Jun2017.csv'
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2-filtered_09Aug2018.csv'

DO_INIT_UPDATE <- TRUE
DO_WRITE_RDATA  <- TRUE
DO_WRITE_NETCDF <- FALSE

filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.par_fixed  <- '../output/par_deoptim_OPT1_04Jul2018.rds'
filename.covariance <- paste('../output/par_LHS2_',appen,'_04Jul2018.RData', sep='')

library(adaptMCMC)
library(ncdf4)
library(foreach)
library(doParallel)

##==============================================================================
## Data
##=====

# requires: data_to_assim (above) and filename.data (above)
source('GEOCARB-2014_getData.R')

# possible filtering out of some data points with too-narrow uncertainties in
# co2 (causing overconfidence in model simulations that match those data points
# well)
# set to +65%, - 30% uncertain range around the central estimate
if(co2_uncertainty_cutoff > 0) {
  co2_halfwidth <- 0.5*(data_calib$co2_high - data_calib$co2_low)
  ind_filter <- which(co2_halfwidth < co2_uncertainty_cutoff)
  ind_remove <- NULL
  for (ii in ind_filter) {
    range_original <- data_calib[ii,'co2_high']-data_calib[ii,'co2_low']
    range_updated  <- data_calib[ii,'co2']*0.95
    if (range_updated > range_original) {
      # update to the wider uncertain range if +65/-30% is wider
      data_calib[ii,'co2_high'] <- data_calib[ii,'co2']*1.65
      data_calib[ii,'co2_low']  <- data_calib[ii,'co2']*0.70
    } else {
      # otherwise, remove
      ind_remove <- c(ind_remove, ii)
    }
  }
  data_calib <- data_calib[-ind_filter,]    # removing all of the possibly troublesome points
  ##data_calib <- data_calib[-ind_remove,]    # remove only those the revised range does not help
}

#===============================================================================


# how many of the n_data_total data points will be in each of the n_shard subsamples?
n_data_total <- nrow(data_calib)
n_data_subsample <- floor(n_data_total/.n_shard)

# TODO -- code with the logical gates from above








# store all of the data_calib subsamples - put all remaining in the last one
data_calib_subsamples <- vector('list', .n_shard)
proxy_types <- as.character(unique(data_calib$proxy_type))
n_proxy <- length(proxy_types)
n_shard_per_proxy <- .n_shard/n_proxy
if(.n_shard==1) {data_calib_subsamples[[1]] <- data_calib} else {
#  ind_remaining <- 1:n_data_total
#  for (s in 1:(.n_shard-1)) {
#    ind_subsample <- sample(ind_remaining, size=n_data_subsample, replace=FALSE)
#    data_calib_subsamples[[s]] <- data_calib[ind_subsample,]
#    ind_remaining <- ind_remaining[-which(ind_remaining %in% ind_subsample)]
#  }
#  data_calib_subsamples[[.n_shard]] <- data_calib[ind_remaining,]
  for (dtype in 1:n_proxy) {
    just_these_data <- data_calib[which(data_calib$proxy_type==proxy_types[dtype]),]
    age_sorted <- sort(just_these_data$age)
    d_age <- floor(length(age_sorted)/n_shard_per_proxy)
    breaks <- age_sorted[d_age*seq(1,(n_shard_per_proxy-1))]
    breaks <- c(0, breaks, 1e4)
    for (sp in 1:n_shard_per_proxy) {
      s <- n_shard_per_proxy*(dtype-1) + sp
      ind_subsample <- which((just_these_data$age >= breaks[sp]) & (just_these_data$age < breaks[sp+1]))
      data_calib_subsamples[[s]] <- just_these_data[ind_subsample,]
    }
  }
}
