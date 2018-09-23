##==============================================================================
## GEOCARB_fit_timevarying_matrix.R
##
## Read centers and standard deviations for the time-varying parameters (vary
## each time step), and put a correlated prior on each time series.
##
## Requires:
## Yields:
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

######
## just sample uncertainty uniformly for the entire time series like:
q_samples <- runif(n=10000, min=0, max=1)
samples_cdf <- qnorm(p=q_samples, mean=par_time_center$Sr[1], sd=sqrt(par_time_covar$Sr[1,1]))
######


#install.packages('mlegp')
library(mlegp)


appen <- 'all'
filename.calibinput <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
source('GEOCARB-2014_parameterSetup.R')

#input

#time_arrays


# par_time0  is the means for each time series (columns) as fcn of time (rows)
# ind_time_err gives the column #s for each of the error estimates for par_time0

ts_names <- names(par_time0)

par_time_center <- vector('list', ncol(par_time0))
names(par_time_center) <- ts_names
for (ts in ts_names) {
  #par_time_center[[ts]] <- diag(par_time0[,ts])
  par_time_center[[ts]] <- par_time0[,ts]
}

time_errors <- time_arrays[,ind_time_err]
par_time_covar <- vector('list', ncol(par_time0))
names(par_time_covar) <- ts_names
for (ts in ts_names) {
  par_time_covar[[ts]] <- diag(time_errors[,paste('e',ts,sep='')]^2)
}

##==============================================================================
## Multivariate normal sampling
##    do a precalibration to get mean and covariance for each time series
##==============================================================================

par_time_center_updated <- par_time_center
par_time_covar_updated <- par_time_covar

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# Distribution fit to each data point in the processing step
#dist <- 'ga'
#dist <- 'be'
#dist <- 'ln'
dist <- 'sn'
source('GEOCARB_fit_likelihood_surface.R')




ts <- ts_names[1]

n_sample <- 10000

X <- mvrnorm(n=n_sample, mu=par_time_center[[ts]], Sigma=par_time_covar[[ts]])





ll_response <- rep(NA, n_sample)

for (ii in 1:n_sample) {

  par_calib_new <- par_calib0
  ind_replace <- which(parnames_calib==ts)
  par_calib_new[ind_replace] <- X[ii,]

  model_out <- model_forMCMC(par_calib=par_calib_new,
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

  ll_response[ii] <- loglikelihood_smoothed(model_out, likelihood_fit, idx_data)
}

ind_drop <- which(is.infinite(ll_response))
X <- X[-ind_drop,]
ll_response <- ll_response[-ind_drop]


par_time_center_updated[[ts]] <- apply(X=X, FUN=mean, MARGIN=2)
par_time_covar_updated[[ts]] <- cov(X)


##fitgp_test <- mlegp(X, ll_response)



##==============================================================================
## Wishart sampling?
##    of the covariance matrix
##==============================================================================

# Good ol' Wikipedia:
# Choice of parameters
#   The least informative, proper Wishart prior is obtained by setting n = p.
#   The prior mean of Wp(V, n) is nV, suggesting that a reasonable choice for V
#    would be (1/n)Σ_0, where Σ_0 is some prior guess for the covariance matrix.

covar_test <- par_time_covar[[1]][1:3,1:3]
#S_wish_test <- solve(covar_test) # scale matrix ("center")
S_wish_test <- (1/nrow(covar_test))*covar_test
v_wish_test <- nrow(S_wish_test) # degrees of freedom

# might need to scale by degrees of freedom...?
samples <- rWishart(n=10000, df=v_wish_test, Sigma=S_wish_test)




##==============================================================================
## End
##==============================================================================
