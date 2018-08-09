##==============================================================================
## sobolTony.R
##
## The `sensitivity` package Sobol routines all choke with NAs, which occur with
## some regularity due to the `failed_runs` in GEOCARB.  So, we just do the
## calculation ourselves for the first-order and total sensitivity indices.
##
## Assumes the parameters are already appropriately scale for use in the model.
##
## Contains:
##  sobolTony -- serial version, calculates indices according to Sobol' 1991
##               verified calculation in `TEST_sobol.R` get results within <1%
##               of the `sensitivity` package results
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


source("sobol_model.R")

#install.packages('foreach')
#install.packages('doParallel')
library(foreach)
library(doParallel)

##==============================================================================
sobolTony <- function(parameters_sampleA, parameters_sampleB, sens,
                      par_fixed, parnames_calib, parnames_fixed,
                      age, ageN, ind_const_calib, ind_time_calib, ind_const_fixed,
                      ind_time_fixed, ind_expected_time, ind_expected_const,
                      iteration_threshold, input, model_ref=NULL, data_calib=NULL,
                      parallel=FALSE, n_core=1, export_names=NULL,
                      n_boot=0, conf=0.9, second=FALSE){

  if(parallel) {sobol_func <- sobol_model_par} else {sobol_func <- sobol_model}

  n_simulations <- nrow(parameters_sampleA)
  p <- ncol(parameters_sampleA)

  # initialize output
  S <- rep(NA, p)
  T <- rep(NA, p)
  CI_S <- mat.or.vec(p, 2)
  CI_T <- mat.or.vec(p, 2)
  colnames(CI_S) <- c(0.5*(1-conf)*100, (conf+0.5*(1-conf))*100)
  colnames(CI_T) <- c(0.5*(1-conf)*100, (conf+0.5*(1-conf))*100)
  if (second) {
    n_s2 <- p*(p-1)*0.5
    S2 <- rep(NA, n_s2)
    CI_S2 <- mat.or.vec(n_s2, 2)
  }

  conf_lo <- 100*0.5*(1-conf)
  conf_hi <- 100-conf_lo

  if (nrow(parameters_sampleB) != n_simulations) {
    stop("samples A and B must have same number of parameters")
  }


  ## ==========================
  ## ESTIMATE MEAN AND VARIANCE
  ## ==========================

  print('Starting estimation of full model mean and variance...')

  # run the model ensemble under sample A
  mA <- sobol_func(parameters_sampleA, sens, par_fixed, parnames_calib,
                   parnames_fixed, age, ageN, ind_const_calib, ind_time_calib,
                   ind_const_fixed, ind_time_fixed, ind_expected_time,
                   ind_expected_const, iteration_threshold, input,
                   model_ref, data_calib, n_core, export_names)

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
  V0 <- mean(mA^2) - m0^2

  print('... Done estimating full model mean and variance.')



  ## ========================================
  ## ESTIMATE FIRST-ORDER SENSITIVITY INDICES
  ## ========================================

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
    mA_i <- mA  # save these as a local temporary copy, so failed runs can be
                # dropped from all and the calculation for Vi can be consistent
                # in # samples, for bootstrap CI
    m_BA <- sobol_func(p_BA, sens, par_fixed, parnames_calib, parnames_fixed,
                       age, ageN, ind_const_calib, ind_time_calib,
                       ind_const_fixed, ind_time_fixed, ind_expected_time,
                       ind_expected_const, iteration_threshold, input,
                       model_ref, data_calib, n_core, export_names)

    # account for possible NA, and extreme values
    idrop <- which(is.na(m_BA) | is.nan(m_BA) | is.infinite(m_BA) | m_BA < -10)
    if (length(idrop) > 0) {
      idrop_all <- c(idrop_all, idrop)
      m_BA <- m_BA[-idrop]
      mA_i <- mA_i[-idrop]
    }
    Vi[i] <- mean(mA_i * m_BA) - m0^2

    # first-order indices are Si = Vi/V
    S[i] <- Vi[i]/V0

    # confidence interval for S[i]
    if (n_boot > 0) {
      # pick n_CI = length(m_BA) runs from m_BA, mA_i, and get
      # Vi_CI = mean(mA_i[ind_CI]*m_BA[ind_CI]) - mean(mA_i[ind_CI])^2
      S_CI <- rep(NA, n_boot)
      for (j in 1:n_boot) {
        ind_CI <- sample(1:length(m_BA), size=length(m_BA), replace=TRUE)
        Vi_CI <- mean(mA_i[ind_CI]*m_BA[ind_CI]) - mean(mA_i[ind_CI])^2
        V0_CI <- mean(mA_i[ind_CI]^2) - mean(mA_i[ind_CI])^2
        S_CI[j] <- Vi_CI/V0_CI
      }
      CI_S[i,] <- quantile(x=S_CI, probs=c(0.01*conf_lo, 0.01*conf_hi), na.rm=TRUE)
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  print('... Done estimating first-order sensitivity indices.')



  ## ===================================================
  ## ESTIMATE SECOND-ORDER SENSITIVITY INDICES ... MAYBE
  ## ===================================================

  if(second) {
    print('Starting estimation of second-order sensitivity indices...')

    # Vik = variance by including parameters i and k
    #    = [(1/N) sum_{j=1}^N m(xjA)*m(x-ikjB , xikjA)] - m0^2
    # xjA   = jth parameter set of sample A
    # x-ikjB = jth parameter set of B, but with i and k replaced by jth value from A
    Vik <- S2 <- rep(NA, n_s2)
    pb <- txtProgressBar(min=0,max=n_s2,initial=0,style=3)
    cnt <- 1
    for (i in 1:(p-1)) {
      for (k in (i+1):p) {
        # replace the ith column of B with that of A
        p_BA <- parameters_sampleB
        p_BA[,c(i,k)] <- parameters_sampleA[,c(i,k)]
        mA_i <- mA
        m_BA <- sobol_func(p_BA, sens, par_fixed, parnames_calib, parnames_fixed,
                           age, ageN, ind_const_calib, ind_time_calib,
                           ind_const_fixed, ind_time_fixed, ind_expected_time,
                           ind_expected_const, iteration_threshold, input,
                           model_ref, data_calib, n_core, export_names)

        # account for possible NA, and extreme values
        idrop <- which(is.na(m_BA) | is.nan(m_BA) | is.infinite(m_BA) | m_BA < -10)
        if (length(idrop) > 0) {
          idrop_all <- c(idrop_all, idrop)
          m_BA <- m_BA[-idrop]
          mA_i <- mA_i[-idrop]
        }

        Vik[cnt] <- mean(mA_i * m_BA) - m0^2

        # second-order indices are Sik = Vik/V - Si - Sk
        S2[cnt] <- Vik[cnt]/V0 - S[i] - S[k]

        # confidence interval for S2[cnt]
        if (n_boot > 0) {
          # pick n_CI = length(m_BA) runs from m_BA, mA_i, and get
          # Vik_CI = mean(mA_i[ind_CI]*m_BA[ind_CI]) - mean(mA_i[ind_CI])^2
          S2_CI <- rep(NA, n_boot)
          for (j in 1:n_boot) {
            ind_CI <- sample(1:length(m_BA), size=length(m_BA), replace=TRUE)
            Vik_CI <- mean(mA_i[ind_CI]*m_BA[ind_CI]) - mean(mA_i[ind_CI])^2
            V0_CI <- mean(mA_i[ind_CI]^2) - mean(mA_i[ind_CI])^2
            # we are not sampling for Si or Sk above so only consider the
            # variation in Vik for the bootstrapping
            S2_CI[j] <- Vik_CI/V0_CI - S[i] - S[k]
          }
          CI_S2[cnt,] <- quantile(x=S2_CI, probs=c(0.01*conf_lo, 0.01*conf_hi), na.rm=TRUE)
        }

        setTxtProgressBar(pb, cnt)
        cnt <- cnt+1
      }
    }
    close(pb)

    print('... Done estimating second-order sensitivity indices.')
  }



  ## ==================================
  ## ESTIMATE TOTAL SENSITIVITY INDICES
  ## ==================================

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
    m_AB <- sobol_func(p_AB, sens, par_fixed, parnames_calib, parnames_fixed,
                       age, ageN, ind_const_calib, ind_time_calib,
                       ind_const_fixed, ind_time_fixed, ind_expected_time,
                       ind_expected_const, iteration_threshold, input,
                       model_ref, data_calib, n_core, export_names)

    # account for possible NA, and extreme values
    idrop <- which(is.na(m_AB) | is.nan(m_AB) | is.infinite(m_AB) | m_AB < -10)
    if (length(idrop) > 0) {
      idrop_all <- c(idrop_all, idrop)
      m_AB <- m_AB[-idrop]
      mA_i <- mA_i[-idrop]
    }
    V_i[i] <- mean(mA_i * m_AB) - m0^2

    # total indices are Ti = 1 - V_i/V
    T[i] <- 1 - (V_i[i]/V0)

    # confidence interval for S[i]
    if (n_boot > 0) {
      # pick n_CI = length(m_AB) runs from m_AB, mA_i, and get
      # V_i_CI = mean(mA_i[ind_CI]*m_AB[ind_CI]) - mean(mA_i[ind_CI])^2
      T_CI <- rep(NA, n_boot)
      for (j in 1:n_boot) {
        ind_CI <- sample(1:length(m_AB), size=length(m_AB), replace=TRUE)
        V_i_CI <- mean(mA_i[ind_CI]*m_AB[ind_CI]) - mean(mA_i[ind_CI])^2
        V0_CI <- mean(mA_i[ind_CI]^2) - mean(mA_i[ind_CI])^2
        T_CI[j] <- 1 - V_i_CI/V0_CI
      }
      CI_T[i,] <- quantile(x=T_CI, probs=c(0.01*conf_lo, 0.01*conf_hi), na.rm=TRUE)
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  print('... Done estimating total sensitivity indices.')


  ## ======
  ## OUTPUT
  ## ======

  out.S <- cbind(S, CI_S)
  colnames(out.S) <- c('S',paste('S.0',conf_lo, sep=''),paste('S.',conf_hi, sep=''))
  rownames(out.S) <- parnames_calib
  out.T <- cbind(T, CI_T)
  colnames(out.T) <- c('T',paste('T.0',conf_lo, sep=''),paste('T.',conf_hi, sep=''))
  rownames(out.T) <- parnames_calib
  if(second) {
    out.S2 <- cbind(S2, CI_S2)
    colnames(out.S2) <- c('S2',paste('S2.0',conf_lo, sep=''),paste('S2.',conf_hi, sep=''))
    S2.names <- NULL
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        S2.names <- rbind(S2.names, c(parnames_calib[i],parnames_calib[j]))
      }
    }
    #out.S2 <- cbind(S2.names, out.S2)
    out <- list(out.S, out.S2, out.T, S2.names, parameters_sampleA, parameters_sampleB, unique(idrop_all))
    names(out) <- c("S", "S2", "T", "S2.names", "pA", "pB", "idrop")
  } else {
    out <- list(out.S, out.T, parameters_sampleA, parameters_sampleB, unique(idrop_all))
    names(out) <- c("S", "T", "pA", "pB", "idrop")
  }

  return(out)
}

##==============================================================================
## End
##==============================================================================
