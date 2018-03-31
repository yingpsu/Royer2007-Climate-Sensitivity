s.out <- readRDS('../output/sobol_alpha0_30Mar2018.rds')

# sample a bunch of parameters
n_sample <- 10000
x1 <- s.out$X1[sample(1:nrow(s.out$X1), n_sample, replace=FALSE),]

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


# x1: as-is
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
sens1 <- model_present1 - mean(model_present1[is.finite(model_present1)])

# x2: T most sensitive parameters as in x1, other 56-T are held at defaults
T <- 5 <- order(s.out$T[,1], decreasing=TRUE)[1:T]
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
sens2 <- model_present2 - mean(model_present2[is.finite(model_present2)])

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
sens3 <- model_present3 - mean(model_present3[is.finite(model_present3)])

# get rid of the bad runs
irem <- which(is.infinite(sens1) | is.infinite(sens2) | is.infinite(sens3))
sens1 <- sens1[-irem]
sens2 <- sens2[-irem]
sens3 <- sens3[-irem]

# correlations
cc <- cor(cbind(sens1,sens2,sens3))

##
##
##
