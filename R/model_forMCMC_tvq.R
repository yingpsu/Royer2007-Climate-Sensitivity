##==============================================================================
## model_forMCMC_tvq.R
##
## Input:
##  par_calib        vector of input parameters. first N_const_calib are the
##                   time-constant parameters. next ageN are the first
##  par_fixed        vector structured like par_calib, but for fixed arrays
##                   (not calibrated)
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
##  parnames_time    time-varying parameter names
##  do_sample_tvq    (logical) turn CDF samples from par_calib into time-varying
##                   parameter arrays?
##  par_time_center  centers of the time varying parameter arrays
##  par_time_stdev   standard deviations of the time varying parameter arrays
##
## Output:
##  age              age, in millinos of years
##  co2              CO2 concentration, ppmv
##  o2               O2 concentration, ppmv
##==============================================================================

model_forMCMC <- function(par_calib, par_fixed, parnames_calib, parnames_fixed, parnames_time,
                          age, ageN, ind_const_calib, ind_time_calib,
                          ind_const_fixed, ind_time_fixed,
                          ind_expected_time, ind_expected_const,
                          iteration_threshold,
                          do_sample_tvq, par_time_center, par_time_stdev) {

  # this takes in two parameter arrays: one that is all of the parameters
  # actually being calibrated, and the other that is the fixed (non-calib)
  # parameters. here, we paste the two of them together for input to the model.

  N_const_total <- length(c(ind_const_calib, ind_const_fixed))
  N_time_total <- length(c(ind_time_calib, ind_time_fixed))/ageN

  # set up the time-constant parameter matrices
  # first length(ind_const_calib) values are the calibration parameters
  # then the fixed values come at the end
  Matrix_56_unordered <- matrix(c(par_calib[ind_const_calib], par_fixed[ind_const_fixed]),
                      nrow=N_const_total, ncol=1)
  rownames(Matrix_56_unordered) <- c( parnames_calib[ind_const_calib],
                            parnames_fixed[ind_const_fixed] )

  # set up the time-varying parameter matrices
  # rows = time, col = different time series


  # HERE is where mapping from CDF probabilities to quantiles should occur

  # initailize
  Matrix_12_unordered <- as.matrix(par_time_center)

  # sampling
  if(do_sample_tvq) {
    for (ts in parnames_time) {
      cdf_val <- par_calib[match(ts, parnames_calib)]
      Matrix_12_unordered[,ts] <- qnorm(p=cdf_val, mean=par_time_center[[ts]], sd=par_time_stdev[[ts]])
    }
  }

#  geoRes <- GEOCARBSULFvolc_forMCMC(Matrix_56_unordered, Matrix_12_unordered, age, ageN)
  geoRes <- run_geocarbF(Matrix_56=Matrix_56_unordered,
                         Matrix_12=Matrix_12_unordered,
                         age=age,
                         ageN=ageN,
                         iteration_threshold=iteration_threshold,
                         ind_expected_time=ind_expected_time,
                         ind_expected_const=ind_expected_const)

  # for r version
  ###GEOCARB_output <- cbind(age, geoRes$CO2, geoRes$O2)
  # for fortran version
  GEOCARB_output <- cbind(geoRes$age, geoRes$CO2_out, geoRes$O2_out)
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
