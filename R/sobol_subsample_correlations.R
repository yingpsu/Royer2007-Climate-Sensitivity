##==============================================================================
## sobol_subsample_correlations.R
##
## (1) Y. Tang et al 2007 (DOI: 10.5194/hess-11-793-2007)
## (2) T.H. Andres 1997 (DOI: 10.1080/00949659708811804)
## (3) and Nossent et al 2011 (application)(DOI: 10.1016/j.envsoft.2011.08.010)
##
## Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

s.out <- readRDS('../output/sobol_alpha0_30Mar2018.rds')

# how many times to do this experiment, and then average/median the results?
n_iter <- 10
n_sample <- 50000

# get the default parameter values
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# get the reference simulation
model_ref <- model_forMCMC(par_calib=par_calib0,
              par_fixed=par_fixed0,
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
              iteration_threshold=iteration_threshold)[,'co2']

#T_test <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55)
T_test <- c(10, 30, 50)
nT <- length(T_test)

corr_s12 <- mat.or.vec(n_iter, nT)
corr_s13 <- mat.or.vec(n_iter, nT)

for (iter in 1:n_iter) {

  # x1: as-is
  # sample a bunch of parameters
  x1 <- s.out$X1[sample(1:nrow(s.out$X1), n_sample, replace=FALSE),]

  model1 <- sapply(1:n_sample, function(ss) {
        model_forMCMC(par_calib=x1[ss,],
                      par_fixed=par_fixed0,
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
  model_present1 <- apply(X=abs(model1-model_ref), MARGIN=2, FUN=sum)
  sens1_default <- model_present1 - mean(model_present1[is.finite(model_present1)])

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
    x2[,ind_small] <- t(replicate(n_sample, par_calib0[ind_small]))
    model2 <- sapply(1:n_sample, function(ss) {
          model_forMCMC(par_calib=x2[ss,],
                        par_fixed=par_fixed0,
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
    model_present2 <- apply(X=abs(model2-model_ref), MARGIN=2, FUN=sum)
    sens2[[tt]] <- model_present2 - mean(model_present2[is.finite(model_present2)])

    # x3: T most sensitive parameters fixed, other 56-T as in x1
    x3 <- x1
    x3[,ind_large] <- t(replicate(n_sample, par_calib0[ind_large]))
    model3 <- sapply(1:n_sample, function(ss) {
          model_forMCMC(par_calib=x3[ss,],
                        par_fixed=par_fixed0,
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
    model_present3 <- apply(X=abs(model3-model_ref), MARGIN=2, FUN=sum)
    sens3[[tt]] <- model_present3 - mean(model_present3[is.finite(model_present3)])

    # get rid of the bad runs
    irem <- which(is.infinite(sens1[[tt]]) | is.infinite(sens2[[tt]]) | is.infinite(sens3[[tt]]))
    sens1[[tt]] <- sens1[[tt]][-irem]
    sens2[[tt]] <- sens2[[tt]][-irem]
    sens3[[tt]] <- sens3[[tt]][-irem]

    # correlations
    corr[[tt]] <- cor(cbind(sens1[[tt]],sens2[[tt]],sens3[[tt]]))
    corr_s12[iter, tt] <- corr[[tt]][1,2]
    corr_s13[iter, tt] <- corr[[tt]][1,3]

  }
}

save(list=c('x1','sens1','sens2','sens3','corr','T_test','corr_s12','corr_s13'), file='sobol_corr.RData')

##
##
##
