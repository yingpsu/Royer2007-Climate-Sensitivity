##==============================================================================
##
## Codes obtained from Github link included in Lenton et al 2018.
## Last accessed on 23 December 2018.
##==============================================================================

library(gdata)

dat <- read.xls('~/codes/COPSE-COPSE-V2.0/code/forcings/D_haq_inversion_2017.xlsx')
age <- dat[,1]
deg <- dat[,2]

# GEOCARB assumes steady state at 10-Myr increments, so the instantaneous
# value at 10-Myr increments is more consistent than averaging
idx_geocarb <- seq(1,length(age),by=20)
age_geocarb <- age[idx_geocarb]
deg_geocarb <- deg[idx_geocarb]

# trim to the times
ibeg <- match(time_arrays[1,'age'], age_geocarb)
iend <- match(time_arrays[nrow(time_arrays),'age'], age_geocarb)
age_geocarb <- age_geocarb[ibeg:iend]
deg_geocarb <- deg_geocarb[ibeg:iend]

# basically, recreation of Lenton et al 2018 Figure 3A
plot(-age, deg, type='l', col='black')
lines(-age_geocarb, deg_geocarb, col='red', lty=1)
lines(-time_arrays[,'age'], time_arrays[,'fSR'], col='red', lty=2)

# fill in and create a Lenton forcing data set (only difference is fSR)
time_arrays_lenton <- time_arrays
time_arrays_lenton[,'fSR'] <- deg_geocarb

filename.out <- "../input_data/GEOCARB_input_arrays_Lenton.csv"
write.csv(time_arrays_lenton, file=filename.out, )

##==============================================================================
## End
##==============================================================================
