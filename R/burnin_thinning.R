##==============================================================================
## burnin_thinning.R
##
## For control experiment (main text)
##
## Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

library('coda')

setwd('~/codes/GEOCARB/R')

today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_parameters <- paste('../output/processed_mcmc_results_',today,'.RData', sep="")
filename_processing <- paste('../output/processing_',today,'.RData', sep="")

load('../output/geocarb_mcmcoutput_unc_16Jun2019sn.RData')

# get and set up the parameters
USE_LENTON_FSR <- FALSE
USE_ROYER_FSR <- TRUE
filename.calibinput <- "../input_data/GEOCARB_input_summaries_calib_unc.csv"
source('GEOCARB-2014_parameterSetup_tvq.R')

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
#[1] 1.004766

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

lmax <- 3000
cmax <- 0.05
maxlag <- 0

for (m in 1:n_node000) {
    for (p in 1:n_parameters) {
        acf_tmp <- acf(chains_burned[[m]][,p], lag.max=lmax, plot=FALSE)
        idx_low <- which(acf_tmp$acf < cmax)
        while (length(idx_low)==0) {
          lmax <- lmax + 200
          acf_tmp <- acf(chains_burned[[m]][,p], lag.max=lmax, plot=FALSE)
          idx_low <- which(acf_tmp$acf < cmax)
        }
        new <- acf_tmp$lag[idx_low[1]]
        if (maxlag < new) {
            print(paste(m,p,"Updating maxlag to",new))
            maxlag <- new
        }
    }
}

chains_burned_thinned <- chains_burned # initialize
if(n_node000 > 1) {
  for (m in 1:n_node000) {
    chains_burned_thinned[[m]] <- chains_burned[[m]][seq(from=1, to=nrow(chains_burned[[m]]), by=maxlag),]
  }
} else {
  chains_burned_thinned <- chains_burned[seq(from=1, to=nrow(chains_burned), by=maxlag),]
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
save.image(file=filename_processing)

##==============================================================================
## End
##==============================================================================
