##==============================================================================
## GEOCARB_optimization_LHS2.R
##
## only the sensitive parameters vary, others fixed at OPT1 values
## Gives covariance for the sensitive parameter set for MCMC transitions
##
## Questions? Tony Wong (anthony.e.wong@coloradoe.edu)
##==============================================================================


##==============================================================================
## need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')
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
## Now, try a pre-calibration with only the sensitive parameters moving, and
## the others fixed at the DE optim values from above (under calib0).
## Should result in an improved estimate of the covariance for MCMC transitions.

n_sample <- n_sample_lhs
n_parameters <- length(parnames_calib)

## draw parameters by Latin Hypercube (Sample)
parameters_lhs <- randomLHS(n_sample, n_parameters)

## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
for (i in 1:n_const_calib) {
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    par_calib[,i] <- qnorm(p=parameters_lhs[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    par_calib[,i] <- qlnorm(p=parameters_lhs[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
    print('ERROR - unknown prior distribution type')
  }
}

tbeg <- proc.time()

model_out <- sapply(1:n_sample, function(ss) {
                    model_forMCMC(par_calib=par_calib[ss,],
                                  par_fixed=par_fixed1,
                                  parnames_calib=parnames_calib,
                                  parnames_fixed=parnames_fixed,
                                  age=age,
                                  ageN=ageN,
                                  ind_const_calib=ind_const_calib,
                                  ind_time_calib=ind_time_calib,
                                  ind_const_fixed=ind_const_fixed,
                                  ind_time_fixed=ind_time_fixed,
                                  ind_expected_time=ind_expected_time,
                                  ind_expected_const=ind_expected_const,
                                  iteration_threshold=iteration_threshold)[,'co2']})

ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < 100) |
      any(model_out[,ss] > 1e4) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}

tend <- proc.time()

parameters_good <- par_calib[-ibad,]
colnames(parameters_good) <- parnames_calib

# report success rate
print(paste('precalibration took ',(tend[3]-tbeg[3])/60,' minutes', sep=''))
print(paste('with success rate of ',nrow(parameters_good),'/',n_sample, sep=''))

# write output file
new_file <- paste(output.dir,'par_LHS2_',appen,'_',today,'.RData', sep='')
print(paste('writing file ',new_file, sep=''))
save(list=c('parameters_good','par_calib','ibad','model_out'), file=new_file)
##==============================================================================


##==============================================================================
## End
##==============================================================================
