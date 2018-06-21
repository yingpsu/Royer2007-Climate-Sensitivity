##==============================================================================
## sobol.R
##
## The `sensitivity` package Sobol routines all choke with NAs, which occur with
## some regularity due to the `failed_runs` in GEOCARB.  So, we just do the
## calculation ourselves for the first-order and total sensitivity indices.
##
## Assumes the parameters are already appropriately scale for use in the model.
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


source("sobol_model.R")

##==============================================================================
sobol <- function(parameters_sampleA, parameters_sampleB, sens,
                  par_fixed, parnames_calib, parnames_fixed,
                  age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                  ind_time_fixed, ind_expected_time, ind_expected_const,
                  iteration_threshold, input, model_ref=NULL, data_calib=NULL){

  n_simulations <- nrow(parameters_sampleA)
  p <- ncol(parameters_sampleA)

  if (nrow(parameters_sampleB) != n_simulations) {
    stop("samples A and B must have same number of parameters")
  }

  print('Starting estimation of full model mean and variance...')

  # run the model ensemble under sample A
  mA <- sobol_model(parameters_sampleA, sens,
                    par_fixed, parnames_calib, parnames_fixed,
                    age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                    ind_time_fixed, ind_expected_time, ind_expected_const,
                    iteration_threshold, input, model_ref, data_calib)

  # account for possible NA, and extreme values
  idrop <- which(is.na(mA) | is.nan(mA) | is.infinite(mA) | mA < -10)
  idrop_all <- idrop
  if (length(idrop) > 0) {
    mA <- mA[-idrop]
    parameters_sampleA <- parameters_sampleA[-idrop,]
    parameters_sampleB <- parameters_sampleB[-idrop,]
  }

  # estimate m0 and V from sampleA
  m0 <- mean(mA)
  V  <- mean(mA^2) - m0^2

  print('... Done estimating full model mean and variance.')

  print('Starting estimation of first-order sensitivity indices...')

  # Vi = variance by including parameter i
  #    = [(1/N) sum_{j=1}^N m(xjA)*m(x-ijB , xijA)] - m0^2
  # xjA   = jth parameter set of sample A
  # x-ijB = jth parameter set of B, but with i replaced by jth value from A
  Vi <- rep(NA, p)
  pb <- txtProgressBar(min=0,max=p,initial=0,style=3)
  for (i in 1:p) {
    # replace the ith column of B with that of A
    p_BA <- parameters_sampleB
    p_BA[,i] <- parameters_sampleA[,i]
    mA_i <- mA
    m_BA <- sobol_model(p_BA, sens,
                        par_fixed, parnames_calib, parnames_fixed,
                        age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                        ind_time_fixed, ind_expected_time, ind_expected_const,
                        iteration_threshold, input, model_ref, data_calib)

    # account for possible NA, and extreme values
    idrop <- which(is.na(m_BA) | is.nan(m_BA) | is.infinite(m_BA) | m_BA < -10)
    if (length(idrop) > 0) {
      idrop_all <- c(idrop_all, idrop)
      m_BA <- m_BA[-idrop]
      mA_i <- mA_i[-idrop]
    }
    Vi[i] <- mean(mA_i * m_BA) - m0^2
    setTxtProgressBar(pb, i)
  }
  close(pb)

  print('... Done estimating first-order sensitivity indices.')

  print('Starting estimation of total sensitivity indices...')

  # V_i = variance by excluding parameter i
  #    = [(1/N) sum_{j=1}^N m(xjA)*m(x-ijA , xijB)] - m0^2
  # xjA   = jth parameter set of sample A
  # x-ijB = jth parameter set of B, but with i replaced by jth value from A
  V_i <- rep(NA, p)
  pb <- txtProgressBar(min=0,max=p,initial=0,style=3)
  for (i in 1:p) {
    # replace the ith column of A with that of B
    p_AB <- parameters_sampleA
    p_AB[,i] <- parameters_sampleB[,i]
    mA_i <- mA
    m_AB <- sobol_model(p_AB, sens,
                        par_fixed, parnames_calib, parnames_fixed,
                        age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                        ind_time_fixed, ind_expected_time, ind_expected_const,
                        iteration_threshold, input, model_ref, data_calib)

    # account for possible NA, and extreme values
    idrop <- which(is.na(m_AB) | is.nan(m_AB) | is.infinite(m_AB) | m_AB < -10)
    if (length(idrop) > 0) {
      idrop_all <- c(idrop_all, idrop)
      m_AB <- m_AB[-idrop]
      mA_i <- mA_i[-idrop]
    }
    V_i[i] <- mean(mA_i * m_AB) - m0^2
    setTxtProgressBar(pb, i)
  }
  close(pb)

  print('... Done estimating total sensitivity indices.')

  # first-order indices are Si = Vi/V
  S <- Vi/V

  # total indices are Ti = 1 - V_i/V
  T <- 1 - (V_i/V)

  out <- list(S, T, unique(idrop_all))
  names(out) <- c("S", "T", "idrop")

  return(out)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
