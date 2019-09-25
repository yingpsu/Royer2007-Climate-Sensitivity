# latin hypercube with the park and royer 2011 parameters, plus stdev
library(lhs)
n_sample <- 10000
set.seed(2019)
parameters_lhs <- randomLHS(n_sample, length(par_calib0))
par_calib <- parameters_lhs  # initialize

n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
colnames(par_calib) <- parnames_calib
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

model_out <- sapply(1:n_sample, function(ss) {
                    model_forMCMC(par_calib=par_calib[ss,],
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
ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < .lower_bound_co2) |
      any(model_out[,ss] > .upper_bound_co2) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}
parameters_good <- par_calib[-ibad,]


llike_out <- rep(NA, n_sample)
for (ss in 1:n_sample) {
  llike_out[ss] <- loglikelihood_smoothed(model_out[,ss], likelihood_fit, idx_data, stdev=par_calib[ss,7])
}


pp <- 7
plot(density(par_calib[idx_bad,pp]), col="coral", lwd=2, main=parnames_calib[pp], xlab='')
lines(density(par_calib[,pp]), col="steelblue", lwd=2)

##
##
##
