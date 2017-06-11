##==============================================================================
## GEOCARB_lhs_driver_script.R
##
## Read CO2 proxy data. Set which data sets you intend to calibrate using.
## Version for running as script on HPC.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

niter_lhs000 <- 1e3   # number of LHS samples per node
n_node000 <- 1         # number of CPUs to use
#setwd('/home/scrim/axw322/codes/GEOCARB/R')
appen <- ''
#rm(list=ls())


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
source('GEOCARB-2014_parameterSetup.R')
##==============================================================================


##==============================================================================
## Calibration parameter prior distributions
##==========================================

# Get model parameter prior distributions
names <- as.character(input$parameter)
bound_lower <- rep(NA, length(names))
bound_upper <- rep(NA, length(names))

ind_neg_inf <- which(input[,'lower_limit']=='_inf')
bound_lower[ind_neg_inf] <- -Inf
bound_lower[setdiff(1:length(names), ind_neg_inf)] <- as.numeric(as.character(input$lower_limit[setdiff(1:length(names), ind_neg_inf)]))
bound_upper <- input$upper_limit

bounds <- cbind(bound_lower, bound_upper)
rownames(bounds) <- as.character(input$parameter)

# only actually need the calibration parameters' bounds, so reformat the bounds
# array to match the vector of calibration parameters
bounds_calib <- mat.or.vec(nr=length(parnames_calib), nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:length(parnames_calib)) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the calibration MCMC call below.
##==============================================================================


##==============================================================================
## Set up the LHS sample
##======================

# need the physical model
source('model_forMCMC.R')

# need the likelihood function (log_like)
source('GEOCARB-2014_calib_likelihood.R')

library(lhs)
niter_lhs <- niter_lhs000

# set up LHS sample
parameters.lhs0 <- randomLHS(niter_lhs, length(parnames_calib))

# initialize
parameters.lhs <- parameters.lhs0

# resample for the actual parameter distirbutions
for (i in 1:length(parnames_calib)) {
#  parameters.lhs[,i] <- qunif(parameters.lhs0[,i], bounds_calib[i,1], bounds_calib[i,2])
  # instead, draw from the priors
  row_num <- match(parnames_calib[i],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    parameters.lhs[,i] <- qnorm(p=parameters.lhs0[,i], mean=input[row_num,'mean'], sd=(0.5*input[row_num,"two_sigma"]))
    # check if out of bounds
    irem <- which(parameters.lhs[,i] < bounds_calib[i,1] | parameters.lhs[,i] > bounds_calib[i,2])
    parameters.lhs <- parameters.lhs[-irem,]
    parameters.lhs0 <- parameters.lhs0[-irem,]
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    parameters.lhs[,i] <- qlnorm(p=parameters.lhs0[,i], meanlog=log(input[row_num,'mean']), sdlog=log(0.5*input[row_num,"two_sigma"]))
    # check if out of bounds
    irem <- which(parameters.lhs[,i] < bounds_calib[i,1] | parameters.lhs[,i] > bounds_calib[i,2])
    parameters.lhs <- parameters.lhs[-irem,]
    parameters.lhs0 <- parameters.lhs0[-irem,]
  } else {print('ERROR - unknown distribution type')}
}
if(length(parnames_calib) > 1) {colnames(parameters.lhs) <- parnames_calib}

# reset the number of ensemble members, now that some that were out of bounds
# have been removed
if(length(parnames_calib) > 1) {
  niter_lhs <- nrow(parameters.lhs)
} else {
  niter_lhs <- length(parameters.lhs)
}
##==============================================================================


##==============================================================================
## Priors - testing
##=================

# fit beta distributions to the parameters (alpha=shape1, beta=shape2)
beta.names <- c('alpha','beta','lb','ub')
prior.beta <- vector('list', length(beta.names)); names(prior.beta) <- beta.names
for (i in 1:length(beta.names)) {
  prior.beta[[i]] <- rep(NA, length(parnames_calib))
}

# use 'qlims' quantiles for upper and lower bounds (replacing inf/ -inf)

##
## >>> NEED TO MESS AROUND WITH THIS AND SEE WHAT WORKS <<<
qlims <- c(0.0001, 0.9999)
## >>> NEED TO MESS AROUND WITH THIS AND SEE WHAT WORKS <<<
##

for (p in 1:length(parnames_calib)) {
  row_num <- match(parnames_calib[p],input$parameter)
  if(input[row_num, 'distribution_type']=='gaussian') {
    prior.beta$lb[p] <- qnorm(p=qlims[1], mean=input[row_num,'mean'], sd=(0.5*input[row_num,"two_sigma"]))
    prior.beta$ub[p] <- qnorm(p=qlims[2], mean=input[row_num,'mean'], sd=(0.5*input[row_num,"two_sigma"]))
  } else if(input[row_num, 'distribution_type']=='lognormal') {
    prior.beta$lb[p] <- qlnorm(p=qlims[1], meanlog=log(input[row_num,'mean']), sdlog=log(0.5*input[row_num,"two_sigma"]))
    prior.beta$ub[p] <- qlnorm(p=qlims[2], meanlog=log(input[row_num,'mean']), sdlog=log(0.5*input[row_num,"two_sigma"]))
  } else {print('ERROR - unknown distribution type')}
}

# use mean and two_sigma from input_sumaries_table to fit alpha and beta
for (p in 1:length(parnames_calib)) {
  row_num <- match(parnames_calib[p],input$parameter)
  mean_beta <- (input[row_num,'mean'] - prior.beta$lb[p]) / (prior.beta$ub[p] - prior.beta$lb[p])
  sd_beta <- (0.5*input[row_num,'two_sigma']) / (prior.beta$ub[p] - prior.beta$lb[p])
  prior.beta$alpha[p] <- ((1-mean_beta)*(mean_beta^2) / (sd_beta^2)) - mean_beta
  prior.beta$beta[p] <- prior.beta$alpha[p]*((1/mean_beta)-1)
}

# this is how we would generate samples from these priors
#tmp <- qbeta(p=parameters.lhs0, shape1=prior.beta$alpha, shape2=prior.beta$beta)*(prior.beta$ub-prior.beta$lb)+prior.beta$lb

# this is how we would map the parameters back to [0,1] and compute beta prior
# (note that we still have parameters.lhs0, so don't really need to do this)
#tmp <- dbeta(x=(parameters.lhs-prior.beta$lb[p])/(prior.beta$ub[p]-prior.beta$lb[p]), shape1=prior.beta$alpha, shape2=prior.beta$beta)

##==============================================================================


##==============================================================================
## Data
##=====

# Read proxy data. Returns "data_calib_all"
library(sn)
source('GEOCARB-2014_getData.R')

# remove the lowest 10 co2 content data points (all paleosols, it turns out)
# (lowest ~40 are all from paleosols, actually)
ind_co2_sort_all <- order(data_calib_all$co2)
data_calib_all <- data_calib_all[-ind_co2_sort_all[1:20], ]


data_tmp <- cbind( c("paleosols" , FALSE),
                   c("alkenones" , FALSE),
                   c("stomata"   , FALSE),
                   c("boron"     , FALSE),
                   c("liverworts", FALSE) )

nproxy <- rep(0, 5)
for (dd in 1:5) {nproxy[dd] <- length(which(data_calib_all$proxy_type == data_to_assim[1,dd]))}

# initialize
llike_lhs <- mat.or.vec(niter_lhs, ncol(data_tmp))
lprior_lhs <- mat.or.vec(niter_lhs, length(parnames_calib))

# save prior values for beta distribution fits
lprior_beta <- mat.or.vec(niter_lhs, length(parnames_calib))

for (dd in 1:ncol(data_tmp)) {

  data_to_assim <- data_tmp
  data_to_assim[2,dd] <- TRUE
  ind_data    <- which(data_to_assim[2,]==TRUE)
  n_data_sets <- length(ind_data)
  ind_assim   <- vector("list",n_data_sets)
  for (i in 1:n_data_sets) {
    ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
  }

  data_calib <- data_calib_all[unlist(ind_assim),]

  # assumption of steady state in-between model time steps permits figuring out
  # which model time steps each data point should be compared against in advance.
  # doing this each calibration iteration would be outrageous!
  # This assumes the model time step is 10 million years, seq(570,0,by=-10). The
  # model will choke later (in calibration) if this is not consistent with what is
  # set within the actual GEOCARB physical model.
  age_tmp <- seq(570,0,by=-10)
  ttmp <- 10*ceiling(data_calib$age/10)
  ind_mod2obs <- rep(NA,nrow(data_calib))
  for (i in 1:length(ind_mod2obs)){
    ind_mod2obs[i] <- which(age_tmp==ttmp[i])
  }

  ##============================================================================
  ## Run the LHS sample
  ##===================

  # run the ensemble
  pb <- txtProgressBar(min=0,max=niter_lhs,initial=0,style=3)
  for (i in 1:niter_lhs) {
    if(length(parnames_calib) > 1) {
      llike_lhs[i,dd] <- log_like(par=parameters.lhs[i,],
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
                           ind_mod2obs=ind_mod2obs
                           )
      if (dd==1) {
        for (p in 1:length(parnames_calib)) {
          row_num <- match(parnames_calib[p],input$parameter)
          if(input[row_num, 'distribution_type']=='gaussian') {
            lprior_lhs[i,p] <- dnorm(x=parameters.lhs[i,p], mean=input[row_num,'mean'], sd=(0.5*input[row_num,"two_sigma"]), log=TRUE)
          } else if(input[row_num, 'distribution_type']=='lognormal') {
            lprior_lhs[i,p] <- dlnorm(x=parameters.lhs[i,p], meanlog=log(input[row_num,'mean']), sdlog=log(0.5*input[row_num,"two_sigma"]), log=TRUE)
          } else {print('ERROR - unknown distribution type')}
          lprior_beta[i,p] <- dbeta(x=(parameters.lhs[i,p]-prior.beta$lb[p])/(prior.beta$ub[p]-prior.beta$lb[p]), shape1=prior.beta$alpha[p], shape2=prior.beta$beta[p], log=TRUE)
        }
      }
    } else {
      llike_lhs[i,dd] <- log_like(par=parameters.lhs[i],
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
                           ind_mod2obs=ind_mod2obs
                           )
      if (dd==1) {
        row_num <- match(parnames_calib,input$parameter)
        if(input[row_num, 'distribution_type']=='gaussian') {
          lprior_lhs[i] <- dnorm(x=parameters.lhs[i], mean=input[row_num,'mean'], sd=(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else if(input[row_num, 'distribution_type']=='lognormal') {
          lprior_lhs[i] <- dlnorm(x=parameters.lhs[i], meanlog=log(input[row_num,'mean']), sdlog=log(0.5*input[row_num,"two_sigma"]), log=TRUE)
        } else {print('ERROR - unknown distribution type')}
        lprior_beta[i] <- dbeta(x=(parameters.lhs[i]-prior.beta$lb)/(prior.beta$ub-prior.beta$lb), shape1=prior.beta$alpha, shape2=prior.beta$beta, log=TRUE)
      }
    } # end if(length(parnames_calib))
    setTxtProgressBar(pb, i)
  } # end loop over i
  close(pb)
} # end loop over dd
llike_joint <- apply(llike_lhs, 1, sum)

# calculate BIC for each proxy set
bic <- rep(NA, ncol(llike_lhs))
for (dd in 1:ncol(llike_lhs)) {
  bic[dd] <- -2*max(llike_lhs[,dd]) + 1*log(nproxy[dd])
}

# calcualte Schwarz criterion
# define refrence set to be paleosols (wlog)
sc <- -2*(bic-bic[1])

# calcualte Bayes factors, using Schwarz criterion as rough approximation (Kass and Raftery 1995)
bf <- exp(sc)

##==============================================================================
# save
save.image(file=paste('GEOCARB_lhs_',appen,'.RData', sep=''))
##==============================================================================


##==============================================================================
## Plot
##=====

# likelihood functions
par(mfrow=c(3,2))
for (dd in 1:ncol(data_tmp)) {
plot(parameters.lhs, llike_lhs[,dd], xlab='deltaT2X', ylab='log-like', main=data_tmp[1,dd])
  points(parameters.lhs[which.max(llike_lhs[,dd])], max(llike_lhs[,dd]), col='red', pch=16, cex=1.6)
}
plot(parameters.lhs, llike_joint, xlab='deltaT2X', ylab='log-like', main='joint')
  points(parameters.lhs[which.max(llike_joint)], max(llike_joint), col='red', pch=16, cex=1.6)


# "posterior" (likelihood functions times priors)
par(mfrow=c(3,2))
for (dd in 1:ncol(data_tmp)) {
plot(parameters.lhs, llike_lhs[,dd]+lprior_lhs, xlab='deltaT2X', ylab='log-like', main=data_tmp[1,dd])
  points(parameters.lhs[which.max(llike_lhs[,dd])], max(llike_lhs[,dd]), col='red', pch=16, cex=1.6)
}
plot(parameters.lhs, llike_joint+lprior_lhs, xlab='deltaT2X', ylab='log-like', main='joint')
  points(parameters.lhs[which.max(llike_joint+lprior_lhs)], max(llike_joint+lprior_lhs), col='red', pch=16, cex=1.6)


# different weighting schemes
wgt_n <- 1/nproxy; wgt_n <- wgt_n/sum(wgt_n)
llike_joint_wgt_n <- apply(llike_lhs*t(replicate(niter_lhs, wgt_n)), 1, sum)

wgt_bic <- min(bic)/bic; wgt_bic <- wgt_bic/sum(wgt_bic)
llike_joint_wgt_bic <- apply(llike_lhs*t(replicate(niter_lhs, wgt_bic)), 1, sum)


plot(parameters.lhs, llike_joint_wgt_n, xlab='deltaT2X', ylab='log-like', main='joint, weighted by N')
  points(parameters.lhs[which.max(llike_joint_wgt_n)], max(llike_joint_wgt_n), col='red', pch=16, cex=1.6)

plot(parameters.lhs, llike_joint_wgt_bic, xlab='deltaT2X', ylab='log-like', main='joint, weighted by BIC')
  points(parameters.lhs[which.max(llike_joint_wgt_bic)], max(llike_joint_wgt_bic), col='red', pch=16, cex=1.6)


##==============================================================================


##==============================================================================
## Write calibrated parameters output file
##========================================

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:length(parnames_calib)){lmax=max(lmax,nchar(parnames_calib[i]))}

## Name the file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters = paste('geocarb_calibratedParameters_',appen,'_',today,'.nc',sep="")

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('geocarb_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
covjump.var <- ncvar_def('covjump', '', list(dim.parameters,dim.parameters), -999)
outnc <- nc_create(filename.parameters, list(parameters.var,parnames.var, covjump.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames_calib)
ncvar_put(outnc, covjump.var, covjump)
nc_close(outnc)

##==============================================================================


##==============================================================================
## End
##==============================================================================
