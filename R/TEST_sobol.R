##==============================================================================
## TEST_sobol.R
##
## Driver for a small test of the sobol routine against those from Sensitvity
## package.  Uses a model that is essentially linear regression.
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


n_sample <- 1000
n_bootstrap <- -1

##==============================================================================
## set up parameters.  first index is slope, second is intercept

parnames <- c('slope', 'intercept')
lb <- c( 0,  0)
ub <- c(10, 10)
n_parameters <- length(parnames)

##==============================================================================
## grid x to evaluate model on, a model function, and true model

mymodel <- function(parameters, xeval){return(parameters[1]*xeval + parameters[2])}

xeval <- seq(0, 5, by=0.1)
p_true <- c(3, 5)
model_true <- mymodel(p_true, xeval)
model_meas <- sapply(1:length(model_true), function(ix) {model_true[ix]+rnorm(1)})

##==============================================================================
## sample parameters from LHS

library(lhs)
parameters <- randomLHS(n=n_sample*2, k=n_parameters)
for (ip in 1:n_parameters) {parameters[,ip] <- parameters[,ip]*(ub[ip]-lb[ip]) + lb[ip]}
ind_sampleA <- sample(x=1:(2*n_sample), size=n_sample, replace=FALSE)
parameters_sampleA <- parameters[ ind_sampleA,]
parameters_sampleB <- parameters[-ind_sampleA,]
parameters_sampleA_df <- data.frame(parameters_sampleA)
parameters_sampleB_df <- data.frame(parameters_sampleB)

##==============================================================================
## define a sobol model for the sensitivity output
## Assumes input of a data frame/matrix-like object of parameters, and returns a
## vector-like object of model sensitivity output.

sobol_model <- function(parameters, xeval, model_meas, df=FALSE) {
  n_simulations <- nrow(parameters)
  if (df) {
    if (is.null(n_simulations)) {
      model_out <- mymodel(parameters=as.numeric(parameters), xeval=xeval)
      # Nash-Sutcliffe efficiency
      sse <- sum( (model_out - model_meas)^2 )
      sst <- sum( (model_meas - mean(model_meas))^2 )
      model_sens <- 1 - sse/sst
    } else {
      # should be done w matrix mult...
      model_out <- sapply(1:n_simulations, function(ss) {
                          mymodel(parameters=as.numeric(parameters[ss,]), xeval=xeval)})
      # Nash-Sutcliffe efficiency
      model_sens <- rep(NA, n_simulations)
      for (ss in 1:n_simulations) {
        sse <- sum( (model_out[,ss] - model_meas)^2 )
        sst <- sum( (model_meas - mean(model_meas))^2 )
        model_sens[ss] <- 1 - sse/sst
      }
    }
  } else {
    if (is.null(n_simulations)) {
      model_out <- mymodel(parameters=parameters, xeval=xeval)
      # Nash-Sutcliffe efficiency
      sse <- sum( (model_out - model_meas)^2 )
      sst <- sum( (model_meas - mean(model_meas))^2 )
      model_sens <- 1 - sse/sst
    } else {
      # should be done w matrix mult...
      model_out <- sapply(1:n_simulations, function(ss) {
                          mymodel(parameters=parameters[ss,], xeval=xeval)})
      # Nash-Sutcliffe efficiency
      model_sens <- rep(NA, n_simulations)
      for (ss in 1:n_simulations) {
        sse <- sum( (model_out[,ss] - model_meas)^2 )
        sst <- sum( (model_meas - mean(model_meas))^2 )
        model_sens[ss] <- 1 - sse/sst
      }
    }
  }
  return(model_sens)
}

##==============================================================================
## this is all from the `sobol.R` routine that we are testing
## So many variables declared are redundant.

sobolTony <- function(parameters_sampleA, parameters_sampleB, xeval, model_meas){

  n_simulations <- nrow(parameters_sampleA)
  p <- ncol(parameters_sampleA)

  if (nrow(parameters_sampleB) != n_simulations) {
    stop("samples A and B must have same number of parameters")
  }

  print('Starting estimation of full model mean and variance...')

  # run the model ensemble under sample A
  mA <- sobol_model(parameters_sampleA, xeval, model_meas)

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
    m_BA <- sobol_model(p_BA, xeval, model_meas)

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
    m_AB <- sobol_model(p_AB, xeval, model_meas)

    V_i[i] <- mean(mA_i * m_AB) - m0^2
    setTxtProgressBar(pb, i)
  }
  close(pb)

  print('... Done estimating total sensitivity indices.')

  # first-order indices are Si = Vi/V
  S <- Vi/V

  # total indices are Ti = 1 - V_i/V
  T <- 1 - (V_i/V)

  out <- list(S, T)#, unique(idrop_all))
  names(out) <- c("S", "T")#, "idrop")
  return(out)
}

##==============================================================================
## run it

s.tony <- sobolTony(parameters_sampleA, parameters_sampleB, xeval, model_meas)

##==============================================================================
## what do the routines in the `sensitivity` package give?

library(sensitivity)

s.sens <- s.tony

# first order indices using the vanilla Sobol' method
s.out <- sobol(model=sobol_model, parameters_sampleA_df, parameters_sampleB_df,
               nboot=0, xeval=xeval, model_meas=model_meas, df=TRUE)
s.sens$S <- s.out$S$original

# total sensitivity indices
s.out <- sobol2002(model=sobol_model, parameters_sampleA_df, parameters_sampleB_df,
                   nboot=0, xeval=xeval, model_meas=model_meas, df=TRUE)
s.sens$T <- s.out$T$original

##==============================================================================
## comparison

# difference in first-order indices
p_err <- rbind((s.tony$S - s.sens$S)/s.sens$S, (s.tony$T - s.sens$T)/s.sens$T)*100

print('Percent difference between Tonys code and Sensitivity package:')
print(paste('  S:', p_err[1,1], p_err[1,2]))
print(paste('  T:', p_err[2,1], p_err[2,2]))

##==============================================================================
## End
##==============================================================================
