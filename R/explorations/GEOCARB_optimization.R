##==============================================================================
## GEOCARB_optimization.R
##
##
## Questions? Tony Wong (anthony.e.wong@coloradoe.edu)
##==============================================================================

rm(list=ls())

##==============================================================================
## all the set-up should be done starting here...

## which parts of the experiment to do?
DO_OPT1 <- FALSE
DO_LHS2 <- TRUE
DO_OPT2 <- FALSE

## set up the evolutionary optimization for parameter estimates
NITER.DEOPTIM <- 5        # number of iterations for evolutionary optimization
NP.DEOPTIM <- 10*2        # population size for DEoptim (do at least 10*[N parameters])
F.DEOPTIM <- 0.8           # as suggested by Storn et al (2006)
CR.DEOPTIM <- 0.9          # as suggested by Storn et al (2006)

## how many latin hypercube precalibration samples? with only sensitive
## parameters varying
n_sample_lhs <- 10000

## where should output be stored?
output.dir <- '../output/'

## cutoff to filter data that has unrealistically narrow uncertainty range
co2_uncertainty_cutoff <- 20

## what sensitivity correlations output to use?
#appen <- 'sig18+GLAC+LIFE'
appen <- 'sig18'
#appen <- 'all'

filename.calib_all <- '../input_data/GEOCARB_input_summaries_calib_all.csv'
filename.calib_sens <- paste('../input_data/GEOCARB_input_summaries_calib_',appen,'.csv', sep='')
filename.par_deoptim <- paste('../output/par_deoptim_OPT1_04Jul2018.rds', sep='')

## ... and ending here

##==============================================================================
## navigate to proper working directory and set up number of cores for possible
## parallelizations

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/GEOCARB/R')
  .Ncore <- 2
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 15  # use multiple cores to process large data?
}

## for file naming
today <- Sys.Date(); today <- format(today,format="%d%b%Y")

##==============================================================================
## bring in relevant libraries

library(sn)
library(DEoptim)
library(lhs)

##==============================================================================
## run the pipeline

if(DO_OPT1) {source('GEOCARB_optimization_OPT1.R')}

if(DO_LHS2) {source('GEOCARB_optimization_LHS2.R')}

if(DO_OPT2) {source('GEOCARB_optimization_OPT2.R')}

##==============================================================================


##==============================================================================
## End
##==============================================================================
