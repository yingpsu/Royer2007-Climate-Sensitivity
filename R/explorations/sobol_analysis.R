##==============================================================================
## sobol_round1_analysis.R
##
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

library(repr)

s.out <- readRDS('../output/sobol_alpha0_NS-n40K-bs10K_25Jun2018.rds')

load('../output/sobol_corr_100.RData')

corr_s12_avg <- apply(corr_s12, 2, mean)
corr_s13_avg <- apply(corr_s13, 2, mean)

n_iter <- nrow(corr_s12)

options(repr.plot.width=7, repr.plot.height=4)
par(mfrow=c(1,2))
plot(corr_s12[1,], xlab='# parameters', ylab='Correlation', main='X1, X2'); for (i in 2:n_iter) {points(corr_s12[i,])}
plot(corr_s13[1,], xlab='# parameters', ylab='Correlation', main='X1, X3'); for (i in 2:n_iter) {points(corr_s13[i,])}


t13 <- T_test[which(corr_s13_avg < 0.05)]
t12 <- T_test[which(corr_s12_avg > 0.9)]

# Find the lowest number of parameters to include that satisfies both constraints
T_min <- intersect(t12, t13)[1]

# Get the parameter names and other input (priors, etc)
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')

# Create a parameters file for this calibration/sensitivity analysis
# this "all" file has 1s for all calibration parameters (constants only), so we
# need to reset the fixed ones to a 0 (or set all to 0 and calib ones to 1)
calib_sig <- read.csv('../input_data/GEOCARB_input_summaries_calib_all.csv')

sens_total <- s.out$T[,1]

ind_large <- order(sens_total, decreasing=TRUE)[1:T_min]
ind_small <- order(sens_total, decreasing=FALSE)[1:(56-T_min)]


# File for significant parameters set
calib_sig$calib <- 0
for (k in 1:length(ind_large)) {
  row <- which(calib_sig$parameter==parnames_calib[ind_large[k]])
  calib_sig$calib[row] <- 1
}
write.csv(x=calib_sig, file=paste('../input_data/GEOCARB_input_summaries_calib_sig',T_min,'.csv',sep=''), row.names=FALSE)

##==============================================================================



##==============================================================================
## End
##==============================================================================
