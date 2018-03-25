##==============================================================================
## sobol_round1_analysis.R
##
## 4 preliminary Sobol' experiments were run, sampling the parameters for the
## precalibration (Latin hypercube with windowing) from one of:
##  1) 25-75 percentile range from each marginal prior distribution,
##  2) 17-34% range
##  3) 5-95% range
##  4) 2.5-97.5% range
##
##
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')

## Functions in other files
source('sobol_functions.R')

# Set number of parameters being analyzed
n_params <- 56

sobol_exp <- c('a50','a34','a10','a05','a0'); n_experiment <- length(sobol_exp)
sobol_files <- vector('list', n_experiment); names(sobol_files) <- sobol_exp
s1st <- vector('list', n_experiment); names(s1st) <- sobol_exp
s1st1 <- vector('list', n_experiment); names(s1st1) <- sobol_exp
ind_s1_sig <- vector('list', n_experiment); names(ind_s1_sig) <- sobol_exp
ind_st_sig <- vector('list', n_experiment); names(ind_st_sig) <- sobol_exp

# Set Sobol indices file names

sobol_files$a50 <- "../output/geocarb_sobol-1-tot_alpha50_17Jan2018.txt"
sobol_files$a34 <- "../output/geocarb_sobol-1-tot_alpha34_18Jan2018.txt"
sobol_files$a10 <- "../output/geocarb_sobol-1-tot_alpha10_18Jan2018.txt"
sobol_files$a05 <- "../output/geocarb_sobol-1-tot_alpha05_18Jan2018.txt"
sobol_files$a0  <- "../output/geocarb_sobol-1-tot_alpha0_20Mar2018.txt"

## Import data from sensitivity analysis
# First- and total-order indices
for (aa in sobol_exp) {
  s1st[[aa]] <- read.csv(sobol_files[[aa]],
                         sep=' ',
                         header=TRUE,
                         nrows = n_params,
                         as.is=c(TRUE,rep(FALSE,5)))
}
parnames.sobol <- s1st[[1]][,1]

## Determine which indices are statistically significant

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
# 'congtr' method means confidence intervals must be fully on the positive side
# and the central estimate must be above 'sig.cutoff'
for (aa in sobol_exp) {
  s1st1[[aa]] <- stat_sig_s1st(s1st[[aa]]
                               ,method="congtr"
                               ,greater=sig.cutoff
                               ,sigCri='either')
}

# of the total sensitivity idnices that are significant for each experiment,
# what lower significance cutoff can we use such that the same set of parameters
# is significant in all experiments?

#st <- s1st1

for (aa in sobol_exp) {
  ind_s1_sig[[aa]] <- which(s1st1[[aa]]$s1_sig==1)
  ind_st_sig[[aa]] <- which(s1st1[[aa]]$st_sig==1)
  #st[[aa]] <- s1st1[[aa]][ind_st_sig[[aa]]]
}

# testing cutoffs for which total sensitivity indices are considered "significant"
cutoff_test <- seq(0.01, 0.3, by=0.005)
n_test <- length(cutoff_test)

# st_sig_exp has 56 rows (1 for each parameter) and n_test columns (at the tested
# levels of significance (total sensitivity index), is the given parameter significant?)
st_sig_exp <- vector('list', n_experiment); names(st_sig_exp) <- sobol_exp
for (aa in sobol_exp) {st_sig_exp[[aa]] <- mat.or.vec(n_params, n_test)}

for (cc in 1:n_test) {
  for (aa in sobol_exp) {
    st_sig_exp[[aa]][s1st1[[aa]]$ST > cutoff_test[cc],cc] <- 1
  }
}

n_st_sig <- mat.or.vec(n_experiment, n_test)
rownames(n_st_sig) <- sobol_exp
colnames(n_st_sig) <- cutoff_test
for (cc in 1:n_test) {
  for (aa in 1:n_experiment) {
    n_st_sig[aa, cc] <- sum(st_sig_exp[[aa]][,cc])
  }
}

pdiff <- rep(-1, n_test)
for (k in 1:n_test) {
    pdiff[k] <- max(n_st_sig[,k])-min(n_st_sig[,k])
}
ind_agree <- which(pdiff <= 1)

##==============================================================================
## Create parameter calibration files for each of 3 sets of parameters:
## 1) S1 = all 31 parameters with ST > 1% for E50 (or ST > 8% for E95, or...)
## 2) S2 = all 21 parameters with S1 > 1% for E50
## 3) S3 = all 13 parameters with ST > 25% for all experiments

ind_s1 <- which(st_sig_exp$a50[,1]==1)
ind_s2 <- which(s1st1$a50$s1_sig==1)
ind_s3 <- which(st_sig_exp$a05[,ind_agree[1]]==1)

# Get the parameter names and other input (priors, etc)
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')

# Create a parameters file for this calibration/sensitivity analysis
# this "all" file has 1s for all calibration parameters (constants only), so we
# need to reset the fixed ones to a 0 (or set all to 0 and calib ones to 1)
calib_s1 <- calib_s2 <- calib_s3 <- read.csv('../input_data/GEOCARB_input_summaries_calib_all.csv')

# File for S1
calib_s1$calib <- 0
for (k in 1:length(ind_s1)) {
  row <- which(calib_s1$parameter==parnames_calib[ind_s1[k]])
  calib_s1$calib[row] <- 1
}
write.csv(x=calib_s1, file='../input_data/GEOCARB_input_summaries_calib_S1.csv', row.names=FALSE)

# File for S2
calib_s2$calib <- 0
for (k in 1:length(ind_s2)) {
  row <- which(calib_s2$parameter==parnames_calib[ind_s2[k]])
  calib_s2$calib[row] <- 1
}
write.csv(x=calib_s2, file='../input_data/GEOCARB_input_summaries_calib_S2.csv', row.names=FALSE)

# File for S3
calib_s3$calib <- 0
for (k in 1:length(ind_s3)) {
  row <- which(calib_s3$parameter==parnames_calib[ind_s3[k]])
  calib_s3$calib[row] <- 1
}
write.csv(x=calib_s3, file='../input_data/GEOCARB_input_summaries_calib_S3.csv', row.names=FALSE)

##==============================================================================




##==============================================================================
## End
##==============================================================================
