##==============================================================================
## scratch_paleosols
##
## What's up with the paleosols?
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================


##==============================================================================
## Data
##=====

# Read proxy data. Returns "data_calib_all"
library(sn)
source('GEOCARB-2014_getData.R')

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , FALSE),
                        c("stomata"   , FALSE),
                        c("boron"     , FALSE),
                        c("liverworts", FALSE) )

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
##==============================================================================


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
source('GEOCARB-2014_parameterSetup.R')
##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distributions
names <- as.character(input$parameter)
bound_lower <- rep(NA, length(names))
bound_upper <- rep(NA, length(names))

ind_neg_inf <- which(input[,'lower_limit']=='_inf')
bound_lower[ind_neg_inf] <- -Inf
bound_lower[setdiff(1:length(names), ind_neg_inf)] <- as.numeric(as.character(input$lower_limit[setdiff(1:length(names), ind_neg_inf)]))
bound_upper <- input$upper_limit

bounds <- cbind(bound_lower, bound_upper)
rownames(bounds) <- as.character(input$parameter)

# only actually need the calibration parameters' bounds, so reformat the bounds
# array to match the vector of calibration parameters
bounds_calib <- mat.or.vec(nr=length(parnames_calib), nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:length(parnames_calib)) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the calibration MCMC call below.
##==============================================================================


##==============================================================================
## Run a 'control' simulation
##===========================

# need the physical model
source('model_forMCMC.R')

# need the likelihood function and prior distributions
source('GEOCARB-2014_calib_likelihood.R')


model_out <- model_forMCMC(par=par_calib0,
                             par_fixed=par_fixed0,
                             parnames_calib=parnames_calib,
                             parnames_fixed=parnames_fixed,
                             age=age,
                             ageN=ageN,
                             ind_const_calib=ind_const_calib,
                             ind_time_calib=ind_time_calib,
                             ind_const_fixed=ind_const_fixed,
                             ind_time_fixed=ind_time_fixed)

model_stdy <- model_out[ind_mod2obs,'co2']
llike <- sum( sapply(1:length(model_stdy), function(i) dsn(x=model_stdy[i],
                       xi=data_calib$xi_co2[i], omega=data_calib$omega_co2[i],
                       alpha=data_calib$alpha_co2[i], log=TRUE)) )

# likelihood function for each data point
llike_all <- sapply(1:length(model_stdy), function(i) dsn(x=model_stdy[i],
                       xi=data_calib$xi_co2[i], omega=data_calib$omega_co2[i],
                       alpha=data_calib$alpha_co2[i], log=TRUE))

# likelihood function sorted from highest (best match) to lowest (worst match)
ind_llike_sort <- rev(order(llike_all))
llike_all_sort <- llike_all[ind_llike_sort]
llike_sort_sum <- cumsum(llike_all_sort

# some analysis:
> cbind(ind_llike_sort, llike_sort_sum)

> data_calib[499:500,]
      X   age age_old age_young  co2 type co2_low co2_high co2_range distribution proxy_type                                                     reference
499 502 210.9      NA        NA -127 mean     -64     -254         1  skew_normal  paleosols Nordt et al., 2015 (expands and refines Atchley et al., 2013)
500 503 210.7      NA        NA -194 mean     -97     -388         1  skew_normal  paleosols Nordt et al., 2015 (expands and refines Atchley et al., 2013)
       xi_co2    omega_co2 alpha_co2
499 -148.3333 2.980934e-15  77.30957
500 -226.3333 7.942354e-16 -52.49290
> itmp <- which(data_calib$co2 < 0)
> length(itmp)
[1] 2
> itmp
[1] 499 500
> data_calib[ind_llike_sort[500:521],]

> ind_co2_sort <- order(data_calib$co2)
> data_calib$co2[ind_co2_sort[1:10]]
 [1] -194 -127    6   10   10   13   15   38   42   45
> width <- data_calib$co2_high - data_calib$co2_low
> ind_width_sort <- order(width)
> ind_width_sort[1:10]

# leads one to believe that it is (1) overconfident co2_low - co2_high range in
# the last about 10 data points here and (2) the negative co2 values that are
# being fit a little wonky. check the fitting later, but for now just throw them
# out and see if things look okay otherwise.
ind_co2_sort <- order(data_calib$co2)
data_calib
##==============================================================================



##==============================================================================
## End
##==============================================================================
