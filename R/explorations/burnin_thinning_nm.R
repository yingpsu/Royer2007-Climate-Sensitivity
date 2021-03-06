##==============================================================================
## burnin_thinning_nm.R
##
## For the normal kernel suppl. experiment, computes the Gelman and Rubin (1992)
## potential scale reduction factor to diagnose Markov chain convergence, chops
## off for burn-in from the set of parallel chains, computes the autocorrelation
## function (ACF) for each parameter of each chain and determines the lag for
## thinning necessary to maintain an autocorrelation maximum of 0.05.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

library('coda')

setwd('~/codes/GEOCARB/R')

today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_parameters <- paste('../output/processed_mcmc_results_normal_',today,'.RData', sep="")

# load up the control experiment results
appen <- "dFpAUsRlMnm"
datestamp <- "18Sep2019"
load(paste('../output/geocarb_mcmcoutput_',appen,'_',datestamp,'.RData', sep=''))

# get and set up the parameters
USE_LENTON_FSR <- FALSE
USE_DT2019_FSR <- TRUE
filename.calibinput <- "../input_data/GEOCARB_input_summaries_calib_all_stdev.csv"
source('parameterSetup_tvq.R')

##==============================================================================
## Burn in
##========

niter_mcmc <- nrow(amcmc_par1[[1]]$samples)
n_parameters <- ncol(amcmc_par1[[1]]$samples)
n_node000 <- length(amcmc_par1)
niter.test <- 0
gr.test <- rep(0, length(niter.test))

if(n_node000 == 1) {
    # don't do GR stats, just cut off first half of chains
    print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
} else if(n_node000 > 1) {
    # this case is FAR more fun
    # accumulate the names of the soon-to-be mcmc objects
    string.mcmc.list <- 'mcmc1'
    for (m in 2:n_node000) {
        string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
    }
    for (i in 1:length(niter.test)) {
        for (m in 1:n_node000) {
            # convert each of the chains into mcmc object
            eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_par1[[m]]$samples[(niter.test[i]+1):niter_mcmc,])', sep='')))
        }
        eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

        gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
    }
} else {print('error - n_node000 < 1 makes no sense')}

# Convergence check:
#> gr.test
#[1] 1.003042

# hack off first ? iterations for burn in
ifirst <- NA
if(n_node000==1) {
  ifirst <- round(0.5*niter_mcmc)
} else {
  gr.max <- 1.1
  for (i in seq(from=length(niter.test), to=1, by=-1)) {
    if( all(gr.test[i:length(gr.test)] < gr.max) ) {ifirst <- niter.test[i]}
  }
}

chains_burned <- NA
if(n_node000 > 1) {
  chains_burned <- vector('list', n_node000)
  for (m in 1:n_node000) {
    chains_burned[[m]] <- amcmc_par1[[m]]$samples[(ifirst+1):niter_mcmc,]
  }
} else {
  chains_burned <- amcmc_out1$samples[(ifirst+1):niter_mcmc,]
}

##==============================================================================
## Thinning
##=========

source("compute_maxlag.R")
maxlags <- compute_maxlag(chains_burned)

chains_burned_thinned <- chains_burned # initialize
if(n_node000 > 1) {
  for (m in 1:n_node000) {
    chains_burned_thinned[[m]] <- chains_burned[[m]][seq(from=1, to=nrow(chains_burned[[m]]), by=maxlags[m]),]
  }
} else {
  chains_burned_thinned <- chains_burned[seq(from=1, to=nrow(chains_burned), by=maxlags),]
}

if (n_node000 > 1) {
  parameters_posterior <- chains_burned_thinned[[1]]
  for (m in 2:n_node000) {
    parameters_posterior <- rbind(parameters_posterior, chains_burned_thinned[[m]])
  }
} else {
  parameters_posterior <- chains_burned_thinned
}

save(parameters_posterior, file=filename_parameters)

##==============================================================================
## End
##==============================================================================
