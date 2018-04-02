##==============================================================================
## sobol_round1_analysis.R
##
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')

## Functions in other files
source('sobol_functions.R')

# Set number of parameters being analyzed
n_params <- 56

#sobol_exp <- c('a50','a34','a10','a05','a0'); n_experiment <- length(sobol_exp)
sobol_exp <- c('a10','a0'); n_experiment <- length(sobol_exp)
sobol_files <- vector('list', n_experiment); names(sobol_files) <- sobol_exp
sobol_files2 <- vector('list', n_experiment); names(sobol_files2) <- sobol_exp
s1st <- vector('list', n_experiment); names(s1st) <- sobol_exp
s1st1 <- vector('list', n_experiment); names(s1st1) <- sobol_exp
s2 <- vector('list', n_experiment); names(s2) <- sobol_exp
s2_sig1 <- vector('list', n_experiment); names(s2_sig1) <- sobol_exp
s2_table <- vector('list', n_experiment); names(s2_table) <- sobol_exp
s2_conf_low <- vector('list', n_experiment); names(s2_conf_low) <- sobol_exp
s2_conf_high <- vector('list', n_experiment); names(s2_conf_high) <- sobol_exp
ind_s1_sig <- vector('list', n_experiment); names(ind_s1_sig) <- sobol_exp
ind_s2_sig <- vector('list', n_experiment); names(ind_s2_sig) <- sobol_exp
ind_st_sig <- vector('list', n_experiment); names(ind_st_sig) <- sobol_exp

## Set Sobol indices file names
sobol_files$a10 <- "../output/geocarb_sobol-1-tot_alpha10_sensL1_31Mar2018.txt"
sobol_files$a0  <- "../output/geocarb_sobol-1-tot_alpha0_sensL1_31Mar2018.txt"
sobol_files2$a10 <- "../output/geocarb_sobol-2_alpha10_sensL1_31Mar2018.txt"
sobol_files2$a0  <- "../output/geocarb_sobol-2_alpha0_sensL1_31Mar2018.txt"

## Import data from sensitivity analysis
# First- and total-order indices
for (aa in sobol_exp) {
  s1st[[aa]] <- read.csv(sobol_files[[aa]],
                         sep=' ',
                         header=TRUE,
                         nrows = n_params,
                         as.is=c(TRUE,rep(FALSE,5)))
  s2_table[[aa]] <- read.csv(sobol_files2[[aa]],
                             sep=' ',
                             header=TRUE,
                             nrows = n_params*(n_params-1)/2,
                             as.is=c(TRUE,rep(FALSE,4)))
  # Convert second-order to upper-triangular matrix
  s2[[aa]] <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2[[aa]][1:(n_params-1), 2:n_params] = upper.diag(s2_table[[aa]]$S2)
  s2[[aa]] <- as.data.frame(s2[[aa]])
  colnames(s2[[aa]]) <- rownames(s2[[aa]]) <- s1st[[aa]]$Parameter

  # Convert confidence intervals to upper-triangular matrix
  s2_conf_low[[aa]] <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2_conf_high[[aa]] <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
  s2_conf_low[[aa]][1:(n_params-1), 2:n_params] = upper.diag(s2_table[[aa]]$S2_conf_low)
  s2_conf_high[[aa]][1:(n_params-1), 2:n_params] = upper.diag(s2_table[[aa]]$S2_conf_high)

  s2_conf_low[[aa]] <- as.data.frame(s2_conf_low[[aa]])
  s2_conf_high[[aa]] <- as.data.frame(s2_conf_high[[aa]])
  colnames(s2_conf_low[[aa]]) <- rownames(s2_conf_low[[aa]]) <- s1st[[aa]]$Parameter
  colnames(s2_conf_high[[aa]]) <- rownames(s2_conf_high[[aa]]) <- s1st[[aa]]$Parameter

}
parnames.sobol <- s1st[[1]][,1]

## Determine which indices are statistically significant

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
# 'congtr' method means confidence intervals must be fully on the positive side
# and the central estimate must be above 'sig.cutoff'
# S2: using the confidence intervals
for (aa in sobol_exp) {
  s1st1[[aa]] <- stat_sig_s1st(s1st[[aa]]
                               ,method="congtr"
                               ,greater=sig.cutoff
                               ,sigCri='either')
  s2_sig1[[aa]] <- stat_sig_s2(s2[[aa]]
                              ,s2_conf_low[[aa]]
                              ,s2_conf_high[[aa]]
                              ,method='congtr'
                              ,greater=sig.cutoff
                              )
}

# what are all the parameters with significant interaction terms?
for (aa in sobol_exp) {
  ind_s1_sig[[aa]] <- which(s1st1[[aa]]$s1_sig==1)
  ind_st_sig[[aa]] <- which(s1st1[[aa]]$st_sig==1)
  #st[[aa]] <- s1st1[[aa]][ind_st_sig[[aa]]]
  # second-order indices need to be looped over all parameters
  ind_s2_sig[[aa]] <- NULL
  for (p in 1:n_params) {
    ind_s2_sig[[aa]] <- c(ind_s2_sig[[aa]], c(p,which(s2_sig1[[aa]][p,] == 1)))
  }
  ind_s2_sig[[aa]] <- unique(ind_s2_sig[[aa]])
}
##==============================================================================


##==============================================================================
## Create parameter calibration files for T significant parameters

# Get the parameter names and other input (priors, etc)
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')

# Create a parameters file for this calibration/sensitivity analysis
# this "all" file has 1s for all calibration parameters (constants only), so we
# need to reset the fixed ones to a 0 (or set all to 0 and calib ones to 1)
calib_sig <- read.csv('../input_data/GEOCARB_input_summaries_calib_all.csv')

# read correaltion experiment output
load('../output/sobol_corr_n10.RData')

corr_s12_avg <- apply(corr_s12, 2, median)
corr_s13_avg <- apply(corr_s13, 2, median)

T <- 25
ind_large <- order(sens_total, decreasing=TRUE)[1:T]
ind_small <- order(sens_total, decreasing=FALSE)[1:(56-T)]

# File for significant parameters set
calib_sig$calib <- 0
for (k in 1:length(ind_large)) {
  row <- which(calib_sig$parameter==parnames_calib[ind_large[k]])
  calib_sig$calib[row] <- 1
}
write.csv(x=calib_sig, file=paste('../input_data/GEOCARB_input_summaries_calib_sig',T,'.csv',sep=''), row.names=FALSE)
##==============================================================================



##==============================================================================
## End
##==============================================================================
