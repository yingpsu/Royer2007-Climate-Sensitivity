##==============================================================================
##
## Input:
##  par              vector of input parameters. first N_const_calib are the
##                   time-constant parameters. next ageN are the first
##  par_fixed        vector structured like par, but for fixed arrays (not calibrated)
##  age              age, in millions of years
##  ageN             number of time steps
##  ind_const_calib  number of calibration parameters constant in time
##  ind_time_calib   number of calibration parameters varying in time (with ageN
##                   different values)
##  ind_const_fixed  number of fixed parameters constant in time
##  ind_time_fixed   number of fixed parameters varying in time (with ageN
##                   different values)
##  parnames_calib   calibration parameter names
##  parnames_fixed   fixed parameter names
##
## Output:
##  age              age, in millinos of years
##  co2              CO2 concentration, ppmv
##  o2               O2 concentration, ppmv
##==============================================================================

source("GEOCARBSULFvolc_forMCMC.R")

model_forMCMC <- function(par, par_fixed, parnames_calib, parnames_fixed,
                          age, ageN, ind_const_calib, ind_time_calib,
						  ind_const_fixed, ind_time_fixed) {

  # this takes in two parameter arrays: one that is all of the parameters
  # actually being calibrated, and the other that is the fixed (non-calib)
  # parameters. here, we paste the two of them together for input to the model.

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
  colnames(Matrix_12) <- c(parnames_calib[ind_time_calib[seq(1, length(par), by=ageN)]],
                           parnames_time_fixed0_vec[seq(from=1, to=(length(ind_time_fixed)), by=ageN)])

  geoRes <- GEOCARBSULFvolc_forMCMC(Matrix_56, Matrix_12, age, ageN)

  #merge results and export summary file to working directory
  GEOCARB_output <- cbind(age, geoRes$CO2, geoRes$O2)
  colnames(GEOCARB_output) <- c('age','co2','o2')

  return(GEOCARB_output)
}

##==============================================================================
## End
##==============================================================================

if(FALSE) { # don't want to write any output - just return the GEOCARB_output
# to the MCMC
#	write.csv(GEOCARB_output, "GEOCARB_output.csv")

  #####In order to run Yawen's MCMC code "cali_1D_YC.R", I need to first
  #####interpolate the model years such that it has the same output as the data
  #1) read in the output of GEOCARB model
  output <- read.csv("GEOCARB_output.csv")
  #2) Interpolate the model output
  baseline <- read.csv("PhanCO2_input.csv")
  fun <- approxfun(output[,2], output[,4])
  output.approx <- fun(baseline[,3])
  #model.year <- output[,2]
  #index=which (model.year == baseline[,1])
  GEOCARB_output_interp <- matrix(0, nrow=944, ncol=2)
  GEOCARB_output_interp[,1] <- baseline[,3]
  GEOCARB_output_interp[,2] <- output.approx
  write.csv(GEOCARB_output_interp, "GEOCARB_output_interp.csv")
  return(list(y=GEOCARB_output_interp))
}


#CO2 <- read.csv("GEOCARB_output_original.csv")
#age_geocarb <- CO2[,1]
#co2 <- CO2[,4]
#plot(age_geocarb,co2,type='l')
