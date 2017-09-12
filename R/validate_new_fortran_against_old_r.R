##==============================================================================
## check of run_geocarbF.R/run_geocarb.f90 against the old R-only version
## (model_forMCMC.R/GEOCARBSULFvolc_forMCMC.R).
##
## Created originally 7 August 2017 by Tony Wong.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/Royer2007-Climate-Sensitivity/R')

##==============================================================================
## Set things up

# Read proxy data. Returns "data_calib_all"
library(sn)
source('GEOCARB-2014_getData.R')

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

ind_data    <- which(data_to_assim[2,]==TRUE)
n_data_sets <- length(ind_data)
ind_assim   <- vector("list",n_data_sets)
for (i in 1:n_data_sets) {
  ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
}

data_calib <- data_calib_all[unlist(ind_assim),]

# assumption of steady state in-between model time steps permits figuring out
# which model time steps each data point should be compared against in advance.
# doing this each calibration iteration would be outrageous!
# This assumes the model time step is 10 million years, seq(570,0,by=-10). The
# model will choke later (in calibration) if this is not consistent with what is
# set within the actual GEOCARB physical model.
age_tmp <- seq(570,0,by=-10)
ttmp <- 10*ceiling(data_calib$age/10)
ind_mod2obs <- rep(NA,nrow(data_calib))
for (i in 1:length(ind_mod2obs)){
  ind_mod2obs[i] <- which(age_tmp==ttmp[i])
}

# Read parameter information, set up the calibration parameters
source('GEOCARB-2014_parameterSetup.R')

# need the physical models
source('model_forMCMC.R')
source('GEOCARBSULFvolc_forMCMC.R')
source('run_geocarbF.R')

# set up some parameters to run at
par <- par_calib0
par_fixed <- par_fixed0
iteration_threshold <- 10

##==============================================================================
## Get output using old R code

model_R <- GEOCARBSULFvolc_forMCMC(Matrix_56, Matrix_12, age, ageN)


##==============================================================================
## Get output using new R/Fortran code

##
## NOTE -- this stuff will go into a similar wrapper vvvvvvvvvvvvvvvvvvvvvv
##

N_const_total <- length(c(ind_const_calib, ind_const_fixed))
N_time_total <- length(c(ind_time_calib, ind_time_fixed))/ageN

# set up the time-constant parameter matrices
# first length(ind_const_calib) values are the calibration parameters
# then the fixed values come at the end
Matrix_56 <- matrix(c(par[ind_const_calib], par_fixed[ind_const_fixed]),
                    nrow=N_const_total, ncol=1)
rownames(Matrix_56) <- c( parnames_calib[ind_const_calib],
                          parnames_fixed[ind_const_fixed] )

# set up the time-varying parameter matrices
# rows = time, col = different time series
Matrix_12 <- matrix(c(par[ind_time_calib], par_fixed[ind_time_fixed]),
                    nrow=ageN, ncol=N_time_total)
colnames(Matrix_12) <- c(unique(parnames_calib[ind_time_calib]), unique(parnames_fixed[ind_time_fixed]))

##
## NOTE -- this stuff will go into a similar wrapper ^^^^^^^^^^^^^^^^^^^^^^
##

model_F <- run_geocarbF(Matrix_56=Matrix_56,
                         Matrix_12=Matrix_12,
                         age=age,
                         ageN=ageN,
                         iteration_threshold=iteration_threshold,
                         ind_expected_time=ind_expected_time,
                         ind_expected_const=ind_expected_const)

##==============================================================================
## Compare output

print(paste('Mean [fortran]-[r] model CO2 bias (ppmv) is: ',mean(model_F$CO2_out-model_R[,2]),sep=''))
print(paste('Max abs([fortran]-[r]) model CO2 bias (ppmv) is: ',max(abs(model_F$CO2_out-model_R[,2])),sep=''))

print(paste('Mean [fortran]-[r] model O2 bias (%) is: ',mean(model_F$O2_out-model_R[,3]),sep=''))
print(paste('Max abs([fortran]-[r]) model O2 bias (%) is: ',max(abs(model_F$O2_out-model_R[,3])),sep=''))

par(mfrow=c(2,2))
plot(-model_R[,1], model_R[,2], col='black', pch=16, xlab='Age [Mya]', ylab='CO2 [ppmv]')
points(-model_F$age, model_F$CO2_out, col='red', pch=16)
plot(-model_R[,1], model_R[,3], col='black', pch=16, xlab='Age [Mya]', ylab='O2 [%]')
points(-model_F$age, model_F$O2_out, col='red', pch=16)
plot(-model_R[,1], model_F$CO2_out-model_R[,2], col='black', pch=16, xlab='Age [Mya]', ylab='F-R CO2 bias [ppmv]')
plot(-model_R[,1], model_F$O2_out-model_R[,3], col='black', pch=16, xlab='Age [Mya]', ylab='F-R O2 bias [%]')

##==============================================================================
## Compare timing
niter <- 10000

tbeg <- proc.time()
for (n in 1:niter) {
  geoRes <- GEOCARBSULFvolc_forMCMC(Matrix_56, Matrix_12, age, ageN)
}
tend <- proc.time()
dt_R <- tend[3]-tbeg[3]



tbeg <- proc.time()
for (n in 1:niter) {
  model_F <- run_geocarbF(Matrix_56=Matrix_56,
                         Matrix_12=Matrix_12,
                         age=age,
                         ageN=ageN,
                         iteration_threshold=iteration_threshold,
                         ind_expected_time=ind_expected_time,
                         ind_expected_const=ind_expected_const)
}
tend <- proc.time()
dt_F <- tend[3]-tbeg[3]


##==============================================================================
## End
##==============================================================================
