##==============================================================================
## model_setup.R
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

library(adaptMCMC)
library(sn)
library(invgamma)

# If using the Foster et al 2017 data set, which proxy sets to assimilate?
#   (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

# do initialization of parameters & covariance matrix from previous calibration?
DO_PARAM_INIT <- TRUE   # if true, set the two file names below
filename.covarinit <- "../output/covar_init_sn-mix_12Aug2019.rds"
filename.paraminit <- "../output/param_init_sn-mix_12Aug2019.rds"
filename.stdevinit <- "../output/geocarb_mcmcoutput_stdevSpinup_07Sep2019.RData"

if (fSR_choice=="PR2011") {
  USE_LENTON_FSR <- FALSE
  USE_DT2019_FSR <- FALSE
} else if (fSR_choice=="LENTON") {
  USE_LENTON_FSR <- TRUE
  USE_DT2019_FSR <- FALSE
} else if (fSR_choice=="ROYER") {
  USE_LENTON_FSR <- FALSE
  USE_DT2019_FSR <- TRUE
} else {
  print("ERROR: unknown fSR_choice")
}

# create appendix tag and file name for this simulation set
appen <- ""
if (data_choice=="PR2011") {appen <- paste(appen,"dP", sep="")} else if (data_choice=="F2017") {appen <- paste(appen,"dF", sep="")}
if (substr(param_choice,1,6)=="PR2011") {appen <- paste(appen,"pP", sep="")} else if (substr(param_choice,1,3)=="all") {appen <- paste(appen,"pA", sep="")}
if (grepl("stdev", param_choice)) {appen <- paste(appen,"U", sep="")}
if (fSR_choice=="PR2011") {appen <- paste(appen,"sO", sep="")} else if (fSR_choice=="LENTON") {appen <- paste(appen,"sL", sep="")} else if (fSR_choice=="ROYER") {appen <- paste(appen,"sR", sep="")}
if (lhood_choice=="unimodal") {appen <- paste(appen,"lU", sep="")} else if (lhood_choice=="mixture") {appen <- paste(appen,"lM", sep="")}
appen <- paste(appen,dist, sep="")
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.mcmc = paste('../output/geocarb_mcmcoutput_',appen,'_',today,'.RData',sep="")

# upper bound from Royer et al 2014 (should be yielding a failed run anyhow)
# lower bound relaxed in light of additional proxy data
.upper_bound_co2 <- 50000
.lower_bound_co2 <- 0

# sampling the time-varying arrays
if (substr(param_choice,1,6)=="PR2011") {
  DO_SAMPLE_TVQ <- FALSE
} else {
  DO_SAMPLE_TVQ <- TRUE
}

# If initializing the parameters and transition matrix, first set up parameters
# as though all were being calibrated, for getting the component names to
# initialize parameters if calibrating using the reduced PR2011 set.
if (DO_PARAM_INIT) {
  filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_all_stdev.csv'
  source('parameterSetup_tvq.R')
  parnames_calib_all <- parnames_calib
}

# Now, set up the parameters as desired. Possibly the same.
filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',param_choice,'.csv', sep='')
source('parameterSetup_tvq.R')
source('model_forMCMC_tvq.R')

#source('run_geocarbF.R')
source('run_geocarbF_unc.R') # version with extra `stdev` uncertainty statistical parameter
##==============================================================================


##==============================================================================
## Initialization of  parameters and transition matrix
##====================================================

if(DO_PARAM_INIT) {
  # initialize the physics parameters with a previous long calibration
  step_mcmc_all <- readRDS(filename.covarinit)
  par_calib0_all <- readRDS(filename.paraminit) # this will overwrite setting stdev above
  # if not using all 69 parameters
  if (substr(param_choice,1,6)=="PR2011") {
    idx_match <- match(parnames_calib, parnames_calib_all)
    par_calib0 <- par_calib0_all[idx_match]
    step_mcmc <- step_mcmc_all[idx_match, idx_match]
  } else {
    par_calib0 <- par_calib0_all
    step_mcmc <- step_mcmc_all
  }
  # initialize the stdev uncertainty parameter
  load(filename.stdevinit)
  par_calib0[match('stdev',parnames_calib)] <- amcmc_out1$samples[amcmc_out1$n.sample]
  step_mcmc[match('stdev',parnames_calib),match('stdev',parnames_calib)] <- amcmc_out1$cov.jump/10
  # using cov.jump/10 so it doesn't move too far, too fast. will adapt anyway.
  # anything that's still NA becomes a 0 and we'll estimate as we go:
  step_mcmc[which(is.na(step_mcmc))] <- 0
}

##==============================================================================


##==============================================================================
## Data
##=====

source('fit_likelihood_surface.R')
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

rm(list=c('bound_lower','bound_upper','bounds'))

##==============================================================================
## End
##==============================================================================
