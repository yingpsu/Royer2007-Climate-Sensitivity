##==============================================================================
## sensitivity_sobol.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
## Using Sobol' Method (variance based)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


##==============================================================================
## Sobol' wrapper - assumes uniform distributions on parameters
geocarb_sobol_ser <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, data_calib, ind_mod2obs, ind_expected_time,
                          ind_expected_const, iteration_threshold) {
  ll_output <- log_like_sensitivity(par_calib_scaled,
                par_fixed=par_fixed0, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input,
                data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold)
  ll_output_centered <- ll_output - mean(ll_output)
  return(ll_output_centered)
}
##==============================================================================


##==============================================================================
## Sobol' wrapper - assumes uniform distributions on parameters
if(FALSE){
geocarb_sobol_par <- function(par_calib_scaled,
                              par_fixed,
                              parnames_calib,
                              parnames_fixed,
                              age,
                              ageN,
                              ind_const_calib,
                              ind_time_calib,
                              ind_const_fixed,
                              ind_time_fixed,
                              input,
                              data_calib,
                              ind_mod2obs,
                              ind_expected_time,
                              ind_expected_const,
                              iteration_threshold,
                              Ncore,
                              alpha) {
}}
geocarb_sobol_par <- function(par_calib_scaled){

  cores=detectCores()
  #cl <- makeCluster(cores[1]-1) #not to overload your computer
  cl <- makeCluster(Ncore)
  print(paste('Starting cluster with ',Ncore,' cores', sep=''))
  registerDoParallel(cl)

  ## initialize output/input
  n_simulations <- nrow(par_calib_scaled)
  #export_names <- c('log_like_sensitivity', 'log_like', 'model_forMCMC', 'run_geocarbF', 'par_fixed0')
  export_names <- c('log_like_sensitivity', 'log_like', 'model_forMCMC', 'run_geocarbF', 'par_fixed0',
                    'parnames_calib', 'parnames_fixed','age','ageN','ind_const_calib',
                    'ind_time_calib', 'ind_const_fixed','ind_time_fixed','input',
                    'data_calib','ind_mod2obs','ind_expected_time','ind_expected_const',
                    'iteration_threshold','Ncore')

  # initialize
  par_calib <- par_calib_scaled
  output <- rep(0, n_simulations)

  # scale the parameters from [0,1] to their prior distributions
  n_const_calib <- length(ind_const_calib)
  n_simulations <- nrow(par_calib_scaled)
  if (n_const_calib > 0) {
    if(n_simulations > 1) {
      for (i in 1:n_const_calib) {
        row_num <- match(parnames_calib[i],input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          par_calib[,i] <- qnorm(p=par_calib_scaled[,ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          par_calib[,i] <- qlnorm(p=par_calib_scaled[,ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
        } else {
          print('ERROR - unknown prior distribution type')
        }
      }
    } else {
      for (i in 1:n_const_calib) {
        row_num <- match(parnames_calib[i],input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          par_calib[i] <- qnorm(p=par_calib_scaled[ind_const_calib[i]], mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          par_calib[i] <- qlnorm(p=par_calib_scaled[ind_const_calib[i]], meanlog=log(input[row_num,"mean"]), sdlog=log(0.5*input[row_num,"two_sigma"]))
        } else {
          print('ERROR - unknown prior distribution type')
        }
      }
    }
  }

if(FALSE){
  # check that none of the normally-distributed parameters are outside any bounds
  # (passes this check, so skipping it since with millions of parameters could
  #  take a while)
  ind_out <- NULL
  for (i in 1:n_const_calib) {
    row_num <- match(parnames_calib[i],input$parameter)
    if(input[row_num, 'distribution_type']=='gaussian') {
      par_min <- qnorm(p=0.5*alpha, mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      par_max <- qnorm(p=(1-0.5*alpha), mean=input[row_num,"mean"], sd=(0.5*input[row_num,"two_sigma"]))
      if (any(par_calib[,i]<par_min | par_calib[,i]>par_max)) {ind_out <- c(ind_out, i)}
    }
  }
}

  finalOutput <- foreach(ii=1:n_simulations,
                         .combine=c,
                         .packages=c('sn'),
                         .export=export_names,
                         .inorder=FALSE) %dopar% {
    dyn.load("../fortran/run_geocarb.so")
    ll_output <- log_like(par_calib=par_calib[ii,],
                          par_fixed=par_fixed0,
                          parnames_calib=parnames_calib,
                          parnames_fixed=parnames_fixed,
                          age=age,
                          ageN=ageN,
                          ind_const_calib=ind_const_calib,
                          ind_time_calib=ind_time_calib,
                          ind_const_fixed=ind_const_fixed,
                          ind_time_fixed=ind_time_fixed,
                          data_calib=data_calib,
                          ind_mod2obs=ind_mod2obs,
                          ind_expected_time=ind_expected_time,
                          ind_expected_const=ind_expected_const,
                          iteration_threshold=iteration_threshold)
    output[ii] <- ll_output
  }
  stopCluster(cl)
  ll_output_centered <- finalOutput - mean(finalOutput)
  return(ll_output_centered)
}
##==============================================================================


##==============================================================================

# getting about 290000 simulations in 2 minutes (# simulations = (p+2)*n_sample)
# (290000 simulations total <=> n_sample=5000)
n_sample <- .n_sample_sobol
n_bootstrap <- .n_bootstrap_sobol
Ncore <- .Ncore

## Sample parameters (need 2 data frames)
parameters_lhs1 <- randomLHS(n_sample, n_parameters)
parameters_lhs2 <- randomLHS(n_sample, n_parameters)

## Trim so you aren't sampling the extreme cases?
alpha <- 0.5
parameters_lhs1 <- (1-alpha)*parameters_lhs1 + 0.5*alpha
parameters_lhs2 <- (1-alpha)*parameters_lhs2 + 0.5*alpha

## Need data frames as input
#parameters_lhs1 <- data.frame(parameters_lhs1)
#parameters_lhs2 <- data.frame(parameters_lhs2)
#colnames(parameters_lhs1) <- colnames(parameters_lhs2) <- parnames_calib

## Actually run the Sobol'
if(ldoparallel ) {
t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_par,
                         parameters_lhs1,
                         parameters_lhs2,
                         scheme='A',
                         nboot=n_bootstrap))
if(FALSE){
  t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_par,
                           parameters_lhs1,
                           parameters_lhs2,
                           scheme='A',
                           nboot=n_bootstrap,
                           par_fixed=par_fixed0, parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                           ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                           input=input,
                           data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                           ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold, Ncore=.Ncore))
}
} else {
  t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_ser,
                           parameters_lhs1,
                           parameters_lhs2,
                           scheme='A',
                           nboot=n_bootstrap,
                           par_fixed=par_fixed0, parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                           ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                           input=input,
                           data_calib=data_calib, ind_mod2obs=ind_mod2obs,
                           ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold))
}

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('max. sensitivity index=',max_sens_ind,' // max. conf int=',max_conf_int,sep=''))

## Write indices, results we'd need to file


today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames_calib,
                                     s.out$S[,1],
                                     s.out$S[,4],
                                     s.out$S[,5],
                                     s.out$T[,1],
                                     s.out$T[,4],
                                     s.out$T[,5]))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                            'S2_conf_high'), nrow=1)
output2.indices <- s.out$S2[,1]
output2.conf1   <- s.out$S2[,4]
output2.conf2   <- s.out$S2[,5]

# 2nd order index names ordered as: (assuming 39 parameters)
# 1. parnames.sobol[1]-parnames.sobol[2]
# 2. parnames.sobol[1]-parnames.sobol[3]
# 3. parnames.sobol[1]-parnames.sobol[4]
# ... etc ...
# 38. parnames.sobol[1]-parnames.sobol[39] << N=2:39 => p1-p[N]
# 39. parnames.sobol[2]-parnames.sobol[3]
# 40. parnames.sobol[2]-parnames.sobol[4]
# 38+37. parnames.sobol[2]-parnames.sobol[39] << N=3:39 => p2-p[N]
# ... etc ...
names2  <- rownames(s.out$S2)
names2a <- rep(NA, length(names2))
names2b <- rep(NA, length(names2))
cnt <- 1
for (i in seq(from=1, to=(length(parnames_calib)-1), by=1)) {         # i = index of first name
    for (j in seq(from=(i+1), to=(length(parnames_calib)), by=1)) {   # j = index of second name
        names2a[cnt] <- parnames_calib[i]
        names2b[cnt] <- parnames_calib[j]
        cnt <- cnt+1
    }
}

output.2nd <- data.frame(cbind( names2a,
                                names2b,
                                output2.indices,
                                output2.conf1,
                                output2.conf2 ))
write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)




##==============================================================================



##==============================================================================
## End
##==============================================================================
