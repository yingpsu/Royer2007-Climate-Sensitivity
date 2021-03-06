##==============================================================================
## sobol_subsample_correlations.R
##
## (1) Y. Tang et al 2007 (DOI: 10.5194/hess-11-793-2007)
## (2) T.H. Andres 1997 (DOI: 10.1080/00949659708811804)
## (3) and Nossent et al 2011 (application)(DOI: 10.1016/j.envsoft.2011.08.010)
##
## Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

setwd('/home/scrim/axw322/codes/GEOCARB/R')

s.out <- readRDS('../output/sobol_sensNS_precal_24Jun2019.rds')

# how many times to do this experiment, and then average/median the results?
n_iter <- 10
n_sample <- 1000
DO_SAMPLE_TVQ <- TRUE  # sample time series uncertainty by CDF parameters?
USE_LENTON_FSR <- FALSE
USE_DT2019_FSR <- TRUE

# get the default parameter values
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_unc.csv'
source('GEOCARB-2014_parameterSetup_tvq.R')

# get the data
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
#co2_uncertainty_cutoff <- 20
# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )
source('GEOCARB-2014_getData.R')

# get the precalibrated parameter samples
# make sure this is consistent with the precalibration results used for the Sobol' analysis
#filename_in <- '../output/geocarb_precalibration_parameters_alpha0_sensL1_01Apr2018.csv'
#parameters_precal <- read.csv(filename_in)
#n_precal <- nrow(parameters_precal)-1
#bandwidths <- parameters_precal[n_precal+1,]
#parameters_precal <- parameters_precal[-(n_precal+1),]

# need the physical model
source('model_forMCMC_tvq.R')
source('run_geocarbF.R')

# get the reference simulation, if you use L1 or L2 norm as sensitivity measure
model_ref <- model_forMCMC(par_calib=par_calib0,
                           par_fixed=par_fixed0,
                           parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed,
                           parnames_time=parnames_time,
                           age=age,
                           ageN=ageN,
                           ind_const_calib=ind_const_calib,
                           ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed,
                           ind_time_fixed=ind_time_fixed,
                           ind_expected_time=ind_expected_time,
                           ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold,
                           do_sample_tvq=DO_SAMPLE_TVQ,
                           par_time_center=par_time_center,
                           par_time_stdev=par_time_stdev)[,"co2"]

T_test <- seq(from=2, to=55, by=1)
#T_test <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55)
#T_test <- c(10, 30, 50)
nT <- length(T_test)

corr_s12 <- mat.or.vec(n_iter, nT)
corr_s13 <- mat.or.vec(n_iter, nT)

tbeg <- proc.time()

for (iter in 1:n_iter) {

  print(paste('iteration:',iter,'/',n_iter))

  # x1: as-is
  # sample a bunch of parameters
  x1 <- s.out$pA[sample(1:nrow(s.out$pA), n_sample, replace=FALSE),]
  #ind_sample <- sample(1:n_precal, n_sample, replace=FALSE)
  #x1 <- parameters_precal[ind_sample,]

  model1 <- sapply(1:n_sample, function(ss) {
        model_forMCMC(par_calib=x1[ss,],
                      par_fixed=par_fixed0,
                      parnames_calib=parnames_calib,
                      parnames_fixed=parnames_fixed,
                      parnames_time=parnames_time,
                      age=age,
                      ageN=ageN,
                      ind_const_calib=ind_const_calib,
                      ind_time_calib=ind_time_calib,
                      ind_const_fixed=ind_const_fixed,
                      ind_time_fixed=ind_time_fixed,
                      ind_expected_time=ind_expected_time,
                      ind_expected_const=ind_expected_const,
                      iteration_threshold=iteration_threshold,
                      do_sample_tvq=DO_SAMPLE_TVQ,
                      par_time_center=par_time_center,
                      par_time_stdev=par_time_stdev)[,'co2']})

  #model_present1 <- apply(X=abs(model1-model_ref), MARGIN=2, FUN=sum)
  #sens1_default <- model_present1 - mean(model_present1[is.finite(model_present1)])
  # Nash-Sutcliffe efficiency
  sens1_default <- rep(NA, n_sample)
  for (ss in 1:n_sample) {
    model_stdy <- model1[ind_mod2obs, ss]
    icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
    sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
    sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
    sens1_default[ss] <- 1 - sse/sst
  }

  sens1 <- vector('list', nT)
  sens2 <- vector('list', nT)
  sens3 <- vector('list', nT)
  corr <- vector('list', nT)

  for (tt in 1:nT) {

    # x1: as is
    sens1[[tt]] <- sens1_default

    # x2: T most sensitive parameters as in x1, other 56-T are held at defaults
    T <- T_test[tt]
    ind_large <- order(s.out$T[,1], decreasing=TRUE)[1:T]
    ind_small <- order(s.out$T[,1], decreasing=FALSE)[1:(56-T)]
    x2 <- x1
    if (length(ind_small)==1) {
      x2[,ind_small] <- replicate(n_sample, par_calib0[ind_small])
    } else {
      x2[,ind_small] <- t(replicate(n_sample, par_calib0[ind_small]))
    }
    model2 <- sapply(1:n_sample, function(ss) {
          model_forMCMC(par_calib=x2[ss,],
                        par_fixed=par_fixed0,
                        parnames_calib=parnames_calib,
                        parnames_fixed=parnames_fixed,
                        parnames_time=parnames_time,
                        age=age,
                        ageN=ageN,
                        ind_const_calib=ind_const_calib,
                        ind_time_calib=ind_time_calib,
                        ind_const_fixed=ind_const_fixed,
                        ind_time_fixed=ind_time_fixed,
                        ind_expected_time=ind_expected_time,
                        ind_expected_const=ind_expected_const,
                        iteration_threshold=iteration_threshold,
                        do_sample_tvq=DO_SAMPLE_TVQ,
                        par_time_center=par_time_center,
                        par_time_stdev=par_time_stdev)[,'co2']})
    #model_present2 <- apply(X=abs(model2-model_ref), MARGIN=2, FUN=sum)
    #sens2[[tt]] <- model_present2 - mean(model_present2[is.finite(model_present2)])
    # Nash-Sutcliffe efficiency
    sens2[[tt]] <- rep(NA, n_sample)
    for (ss in 1:n_sample) {
      model_stdy <- model2[ind_mod2obs, ss]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
      sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
      sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
      sens2[[tt]][ss] <- 1 - sse/sst
    }

    # x3: T most sensitive parameters fixed, other 56-T as in x1
    x3 <- x1
    x3[,ind_large] <- t(replicate(n_sample, par_calib0[ind_large]))
    model3 <- sapply(1:n_sample, function(ss) {
          model_forMCMC(par_calib=x3[ss,],
                        par_fixed=par_fixed0,
                        parnames_calib=parnames_calib,
                        parnames_fixed=parnames_fixed,
                        parnames_time=parnames_time,
                        age=age,
                        ageN=ageN,
                        ind_const_calib=ind_const_calib,
                        ind_time_calib=ind_time_calib,
                        ind_const_fixed=ind_const_fixed,
                        ind_time_fixed=ind_time_fixed,
                        ind_expected_time=ind_expected_time,
                        ind_expected_const=ind_expected_const,
                        iteration_threshold=iteration_threshold,
                        do_sample_tvq=DO_SAMPLE_TVQ,
                        par_time_center=par_time_center,
                        par_time_stdev=par_time_stdev)[,'co2']})
    #model_present3 <- apply(X=abs(model3-model_ref), MARGIN=2, FUN=sum)
    #sens3[[tt]] <- model_present3 - mean(model_present3[is.finite(model_present3)])
    # Nash-Sutcliffe efficiency
    sens3[[tt]] <- rep(NA, n_sample)
    for (ss in 1:n_sample) {
      model_stdy <- model3[ind_mod2obs, ss]
      icomp <- which(is.finite(model_stdy) & !is.na(model_stdy))  # only compare valid values
      sse <- sum( (model_stdy[icomp] - data_calib$co2[icomp])^2 )
      sst <- sum( (data_calib$co2[icomp] - mean(data_calib$co2[icomp]))^2 )
      sens3[[tt]][ss] <- 1 - sse/sst
    }

    # get rid of the bad runs
    irem <- which(is.infinite(sens1[[tt]]) | is.infinite(sens2[[tt]]) | is.infinite(sens3[[tt]]) |
                  sens1[[tt]] < -10 | sens2[[tt]] < -10 | sens3[[tt]] < -10)
    if (length(irem)>0) {
      sens1[[tt]] <- sens1[[tt]][-irem]
      sens2[[tt]] <- sens2[[tt]][-irem]
      sens3[[tt]] <- sens3[[tt]][-irem]
    }

    # correlations
    corr[[tt]] <- cor(cbind(sens1[[tt]],sens2[[tt]],sens3[[tt]]))
    corr_s12[iter, tt] <- corr[[tt]][1,2]
    corr_s13[iter, tt] <- corr[[tt]][1,3]

  }
}

sens_total <- s.out$T

save(list=c('x1','sens1','sens2','sens3','corr','T_test','corr_s12','corr_s13','sens_total'), file='sobol_corr.RData')

tend <- proc.time()
print(paste(n_iter,' iterations took ',(tend-tbeg)[3]/60,' minutes total', sep=''))

##=============================================================================
## End
##=============================================================================
