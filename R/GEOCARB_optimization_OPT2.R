##==============================================================================
## GEOCARB_optimization_OPT2.R
##
## only the sensitive parameters, others fixed at OPT1 values
## Gives initial parameter estimates for MCMC
##
## Questions? Tony Wong (anthony.e.wong@coloradoe.edu)
##==============================================================================


##==============================================================================
## Data
##=====

source('GEOCARB-2014_getData.R')

# remove the lowest [some number] co2 content data points (all paleosols, it turns out)
# (lowest ~40 are all from paleosols, actually)
#ind_co2_sort_all <- order(data_calib_all$co2)
#n_cutoff <- length(which(data_calib_all$co2 < quantile(data_calib_all$co2, 0.01)))
#data_calib_all <- data_calib_all[-ind_co2_sort_all[1:n_cutoff], ]

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

## Read parameter information, set up the calibration parameters
filename.calibinput <- filename.calib_sens
source('GEOCARB-2014_parameterSetup.R')

## set fixed parameters at optimized values
if(!DO_OPT1) {par_deoptim <- readRDS(filename.par_deoptim)}
par_fixed1 <- par_fixed0
names_replace <- parnames_fixed[ind_const_fixed]
for (name in parnames_fixed[ind_const_fixed]) {
  ind_out <- match(name, parnames_fixed)
  ind_in  <- match(name, names(par_deoptim))
  par_fixed1[ind_out] <- par_deoptim[ind_in]
}
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
## Run the optimization
##=====================

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# need the likelihood function and prior distributions
source('GEOCARB-2014_calib_likelihood.R')

# set up and run the actual optimization
p0.deoptim <- par_calib0             # initialize optimized initial parameters
niter.deoptim <- NITER.DEOPTIM       # number of iterations for evolutionary optimization
NP.deoptim <- NP.DEOPTIM             # population size for DEoptim (do at least 10*[N parameters])
F.deoptim <- F.DEOPTIM               # as suggested by Storn et al (2006)
CR.deoptim <- CR.DEOPTIM             # as suggested by Storn et al (2006)

# get some range for each parameter prior distribution, to use as hard bounds
bound.lower.deoptim <- rep(NA, length(par_calib0))
bound.upper.deoptim <- rep(NA, length(par_calib0))
p_lb <- 0.005
p_ub <- 0.995

for (ip in 1:length(par_calib0)) {
  row_num <- match(parnames_calib[ip],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    bound.lower.deoptim[ip] <- qnorm(p=p_lb, mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
    bound.upper.deoptim[ip] <- qnorm(p=p_ub, mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    bound.lower.deoptim[ip] <- qlnorm(p=p_lb, meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
    bound.upper.deoptim[ip] <- qlnorm(p=p_ub, meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}

tbeg=proc.time()

outDEoptim <- DEoptim(neg_log_post, bound.lower.deoptim, bound.upper.deoptim,
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
        par_fixed=par_fixed1, parnames_calib=parnames_calib,
        parnames_fixed=parnames_fixed, age=age, ageN=ageN,
        ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
        ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
        input=input, time_arrays=time_arrays, bounds_calib=bounds_calib,
        data_calib=data_calib, ind_mod2obs=ind_mod2obs,
        ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
        iteration_threshold=iteration_threshold)

tend=proc.time()
print(paste('Took ',(tend-tbeg)[3]/60,' minutes', sep=''))

par_deoptim = outDEoptim$optim$bestmem
names(par_deoptim) <- parnames_calib

new_file <- paste(output.dir,'par_deoptim_OPT2_',appen,'_',today,'.rds', sep='')
print(paste('writing file ',new_file, sep=''))
saveRDS(par_deoptim, file=new_file)
##==============================================================================


##==============================================================================
## End
##==============================================================================
