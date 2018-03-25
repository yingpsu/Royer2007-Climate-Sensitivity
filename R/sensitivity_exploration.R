##==============================================================================
## exploration_lhs_sensitivity.R
##==============================================================================



## Sample one-at-a-time and examine the impacts on simulation
names_test <- c('LIFE','ACT','GYM','FERT','deltaT2X')
n_test <- length(names_test)
n_sample <- 1000

parameters_lhs <- randomLHS(n_sample, n_test)
colnames(parameters_lhs) <- names_test

## Trim so you aren't sampling the extreme cases?
alpha <- 1-1
parameters_lhs <- (1-alpha)*parameters_lhs + 0.5*alpha

## Fill in par_calib with all defaults, then we'll overwrite the n_test parameters
par_calib <- vector('list', n_test); names(par_calib) <- names_test
model_out <- vector('list', n_test); names(model_out) <- names_test
igood     <- vector('list', n_test); names(igood)     <- names_test

for (ptest in names_test) {

  par_calib[[ptest]] <- t(replicate(n_sample, par_calib0))

  ## scale up to the actual parameter distributions
  row_num <- match(ptest, input$parameter)
  col_num <- match(ptest, parnames_calib)
  if(input[row_num, 'distribution_type']=='gaussian') {
    par_calib[[ptest]][,col_num] <- qnorm(p=parameters_lhs[,ptest], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    par_calib[[ptest]][,col_num] <- qlnorm(p=parameters_lhs[,ptest], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {
      print('ERROR - unknown prior distribution type')
  }

  ## Run the model with all these different combos
  model_out[[ptest]] <- sapply(1:n_sample, function(ss) {
      model_forMCMC(par_calib=par_calib[[ptest]][ss,],
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

  igood$ptest <- NULL
  for (ss in 1:n_sample) {
    if( all(model_out[[ptest]][,ss] > 100) &
        all(model_out[[ptest]][,ss] < 1e5) &
        model_out[[ptest]][58,ss] > 280 & model_out[[ptest]][58,ss] < 400) {
      igood[[ptest]] <- c(igood[[ptest]],ss)
    }
  }
}

## Plot all simulations
par(mfrow=c(3,2))
for (ptest in names_test) {
  plot(-age, model_out[[ptest]][,1], type='l', ylim=c(0,10000))
  for (k in 1:n_sample) {lines(-age, model_out[[ptest]][,k], type='l')}
  for (d in 1:length(data_calib$age)) {points(-data_calib$age[d], data_calib$co2[d], col='red', pch=16)}
  lines(-age, model_ref, col='purple', lw=2)
}

############### OLD #####################

# plot "good" simulations
plot(-age, model_out[,igood[1]], type='l', ylim=c(0,10000))
for (k in igood) {lines(-age, model_out[,k], type='l')}
# add data points
for (d in 1:length(data_calib$age)) {points(-data_calib$age[d], data_calib$co2[d], col='red', pch=16)}
# add reference model
lines(-age, model_ref, col='purple', lw=3)


# plot all simulations
plot(-age, model_out[,1], type='l', ylim=c(0,10000))
for (k in 1:n_sample) {lines(-age, model_out[,k], type='l')}
# add data points
for (d in 1:length(data_calib$age)) {points(-data_calib$age[d], data_calib$co2[d], col='red', pch=16)}
# add reference model
lines(-age, model_ref, col='purple', lw=3)

##==============================================================================
## End
##==============================================================================
