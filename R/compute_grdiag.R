##==============================================================================
## compute_grdiag.R
##
## Compute the Gelman and Rubin (1992) potential scale reduction factor for a
## set of Markov chains.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

compute_grdiag <- function(chains, ibeg) {

  string.mcmc.list <- 'mcmc1'
  for (m in 2:length(chains)) {
      string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
  }
  for (m in 1:length(chains)) {
      # convert each of the chains into mcmc object
      eval(parse(text=paste('mcmc',m,' <- as.mcmc(chains[[m]][(ibeg+1):nrow(chains[[m]]),])', sep='')))
  }
  eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))
  gr.test <- as.numeric(gelman.diag(mcmc_chain_list, autoburnin=FALSE)[2])
  return(gr.test)
}
