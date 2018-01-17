icol_st##==============================================================================
## sensitivity_geocarb.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())

co2_uncertainty_cutoff <- 20

n_sample <- 2e5


if(Sys.info()['nodename']=='Tonys-MacBook-Pro.local') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
  .Ncore <- 2
} else {
  # assume on Napa cluster
  machine <- 'remote'
  setwd('/home/scrim/axw322/codes/GEOCARB/R')
  .Ncore <- 15  # use multiple cores to process large data?
}



##==============================================================================
## Data
##=====

source('GEOCARB-2014_getData.R')

# remove the lowest [some number] co2 content data points (all paleosols, it turns out)
# (lowest ~40 are all from paleosols, actually)
#ind_co2_sort_all <- order(data_calib_all$co2)
#n_cutoff <- length(which(data_calib_all$co2 < quantile(data_calib_all$co2, 0.01)))
#data_calib_all <- data_calib_all[-ind_co2_sort_all[1:n_cutoff], ]

# Which proxy sets to assimilate? (set what you want to "TRUE", others to "FALSE")
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )

ind_data    <- which(data_to_assim[2,]==TRUE)
n_data_sets <- length(ind_data)
ind_assim   <- vector("list",n_data_sets)
for (i in 1:n_data_sets) {
  ind_assim[[i]] <- which(as.character(data_calib_all$proxy_type) == data_to_assim[1,ind_data[i]])
}

data_calib <- data_calib_all[unlist(ind_assim),]

# possible filtering out of some data points with too-narrow uncertainties in
# co2 (causing overconfidence in model simulations that match those data points
# well)
# set to +65%, - 30% uncertain range around the central estimate
if(co2_uncertainty_cutoff > 0) {
  co2_halfwidth <- 0.5*(data_calib$co2_high - data_calib$co2_low)
  ind_filter <- which(co2_halfwidth < co2_uncertainty_cutoff)
  ind_remove <- NULL
  for (ii in ind_filter) {
    range_original <- data_calib[ii,'co2_high']-data_calib[ii,'co2_low']
    range_updated  <- data_calib[ii,'co2']*0.95
    if (range_updated > range_original) {
      # update to the wider uncertain range if +65/-30% is wider
      data_calib[ii,'co2_high'] <- data_calib[ii,'co2']*1.65
      data_calib[ii,'co2_low']  <- data_calib[ii,'co2']*0.70
    } else {
      # otherwise, remove
      ind_remove <- c(ind_remove, ii)
    }
  }
  data_calib <- data_calib[-ind_remove,]
}

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
##==============================================================================


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
source('GEOCARB-2014_parameterSetup.R')
n_parameters <- length(parnames_calib)
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
bounds_calib <- mat.or.vec(nr=n_parameters, nc=2)
colnames(bounds_calib) <- c('lower','upper')
rownames(bounds_calib) <- parnames_calib
for (i in 1:n_parameters) {
  bounds_calib[i,'lower'] <- bounds[parnames_calib[i],'bound_lower']
  bounds_calib[i,'upper'] <- bounds[parnames_calib[i],'bound_upper']
}

# Other characteristics are in 'input' and 'time_arrays', which are fed into
# the sensitivity analysis call below.
##==============================================================================


##==============================================================================
## Set up for optimization
##========================

# need the physical model
source('model_forMCMC.R')
source('run_geocarbF.R')

# need the likelihood function and prior distributions
source('GEOCARB-2014_calib_likelihood.R')

# needed packages
#install.packages('sensitivity')
#install.packages('sn')
#install.packages('lhs')
#install.packages('foreach')
#install.packages('doParallel')
library(sensitivity)
library(sn)
library(lhs)
library(foreach)
library(doParallel)
##==============================================================================


##==============================================================================
## Draw parameters by latin hypercube
##===================================

parameters_lhs <- randomLHS(n_sample, n_parameters)


## Trim so you aren't sampling the extreme cases?
alpha <- .5
parameters_lhs <- (1-alpha)*parameters_lhs + 0.5*alpha


## scale up to the actual parameter distributions
n_const_calib <- length(ind_const_calib)
par_calib <- parameters_lhs  # initialize
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
##==============================================================================


##==============================================================================
## Run precalibration simulations
##===============================

model_out <- sapply(1:n_sample, function(ss) {
    model_forMCMC(par_calib=par_calib[ss,],
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

#source('precalibration_doruns.R')

ibad <- NULL
for (ss in 1:n_sample) {
  if( any(model_out[,ss] < 100) |
      any(model_out[,ss] > 1e5) |
      model_out[58,ss] < 280 | model_out[58,ss] > 400) {
    ibad <- c(ibad,ss)
  }
}

parameters_good <- par_calib[-ibad,]
colnames(parameters_good) <- parnames_calib
##==============================================================================


##==============================================================================
## Fit KDEs to sample from for each parameter
##===========================================

pdf_all <- vector('list', n_parameters)
for (pp in 1:n_parameters){
  tmp <- density(parameters_good[,pp], kernel='gaussian', n=100)
  pdf_all[[pp]] <- tmp; names(pdf_all)[pp] <- parnames_calib[pp]
}

## Write a CSV file with the successful parameter combinations and bandwidths
bandwidths <- rep(NA,n_parameters)
for (pp in 1:n_parameters){
  bandwidths[pp] <- pdf_all[[pp]]$bw
}
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_out <- paste('geocarb_precalibration_parameters_',today,'.csv', sep='')
write.csv(rbind(parameters_good, bandwidths), file=filename_out, row.names=FALSE)
##==============================================================================


##==============================================================================
## Function for sampling from KDEs (inflating a LHS?)
##===================================================

## Read KDE results file, separate into parameters and the bandwidths
filename_in <- 'geocarb_precalibration_parameters_16Jan2018.csv'
parameters_node <- read.csv(filename_in)
n_node <- nrow(parameters_node)-1
bandwidths <- parameters_node[n_node+1,]
parameters_node <- parameters_node[-(n_node+1),]

kde_sample <- function(n_sample, nodes, bandwidths) {
  # preliminaries
  n_node <- nrow(nodes)
  n_par <- length(bandwidths)
  if(n_sample > n_node) {print('ERROR: n_sample > n_node')}

  # choose the node rows out of array `nodes`
  i_sample <- sample(x=1:n_node, size=n_sample, replace=FALSE)

  # sample normal random from around each of the parameter nodes from that row
  # (this achieves joint sampling)
  par_sample <- t(sapply(1:n_sample, function(i) {
       rnorm(n=n_parameters, mean=as.numeric(parameters_node[i,]), sd=as.numeric(bandwidths))}))

  return(par_sample)
}
##==============================================================================



##==============================================================================

## Sobol' wrapper - assumes uniform distributions on parameters
geocarb_sobol_co2_ser <- function(par_calib_scaled, par_fixed, parnames_calib,
                          parnames_fixed, age, ageN, ind_const_calib,
                          ind_time_calib, ind_const_fixed, ind_time_fixed,
                          input, ind_expected_time,
                          ind_expected_const, iteration_threshold) {
  co2_output <- sensitivity_co2(par_calib_scaled, l_scaled=TRUE,
                par_fixed=par_fixed, parnames_calib=parnames_calib,
                parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                input=input, ind_expected_time=ind_expected_time,
                ind_expected_const=ind_expected_const,
                iteration_threshold=iteration_threshold)
#  co2_output_centered <- co2_output - mean(co2_output, na.rm=TRUE)
  co2_output_centered <- co2_output - mean(co2_output[is.finite(co2_output)])
  co2_output_centered[is.infinite(co2_output)] <- 0
  return(co2_output_centered)
}

n_sample <- 500
n_bootstrap <- 5000
Ncore <- 1

source('GEOCARB_sensitivity_co2.R')

## Sample parameters (need 2 data frames)
#parameters_sample1 <- kde_sample(n_sample, parameters_node, bandwidths)
#parameters_sample2 <- kde_sample(n_sample, parameters_node, bandwidths)
## Or take directly from the precalibration?
n_half <- floor(0.5*nrow(parameters_node))
parameters_sample1 <- parameters_node[1:n_half,]
parameters_sample2 <- parameters_node[(n_half+1):(2*n_half),]
colnames(parameters_sample1) <- colnames(parameters_sample2) <- parnames_calib

## Actually run the Sobol'
t.out <- system.time(s.out <- sobolSalt(model=geocarb_sobol_co2_ser,
                           parameters_sample1,
                           parameters_sample2,
                           scheme='A',
                           nboot=n_bootstrap,
                           par_fixed=par_fixed0, parnames_calib=parnames_calib,
                           parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                           ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                           ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                           input=input,
                           ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                           iteration_threshold=iteration_threshold))

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('0.1 * max. sensitivity index=',max_sens_ind,' // max. conf int=',max_conf_int,sep=''))
##==============================================================================



##==============================================================================
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

## load relevant output files

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
plotdir <- '../figures/'

# Set number of parameters being analyzed
n_params <- 56
# Set Sobol indices file names
Sobol_file_1 <- "../output/geocarb_sobol-1-tot_17Jan2018.txt"
Sobol_file_2 <- "../output/geocarb_sobol-2_17Jan2018.txt"
s.out <- readRDS('sobol_16Jan2018.rds')

##==============================================================================




##==============================================================================

## some analysis

## Libraries----
library(RColorBrewer) # good color palettes
library(graphics)     # used when plotting polygons
library(plotrix)      # used when plotting circles

## Functions in other files
source('sobol_functions.R')
source('colorblindPalette.R')

## Import data from sensitivity analysis
# First- and total-order indices
s1st <- read.csv(Sobol_file_1,
                  sep=' ',
                  header=TRUE,
                  nrows = n_params,
                  as.is=c(TRUE,rep(FALSE,5)))

parnames.sobol <- s1st[,1]

# Import second-order indices
s2_table <- read.csv(Sobol_file_2,
               sep=' ',
               header=TRUE,
               nrows = n_params*(n_params-1)/2,
               as.is=c(TRUE,rep(FALSE,4)))

# Convert second-order to upper-triangular matrix
s2 <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
s2 <- as.data.frame(s2)
colnames(s2) <- rownames(s2) <- s1st$Parameter

# Convert confidence intervals to upper-triangular matrix
s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

s2_conf_low <- as.data.frame(s2_conf_low)
s2_conf_high <- as.data.frame(s2_conf_high)
colnames(s2_conf_low) <- rownames(s2_conf_low) <- s1st$Parameter
colnames(s2_conf_high) <- rownames(s2_conf_high) <- s1st$Parameter

####################################
# Determine which indices are statistically significant

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
s1st1 <- stat_sig_s1st(s1st
                      ,method="congtr"
                      ,greater=sig.cutoff
                      ,sigCri='either')

# S1 & ST: using greater than a given value
#s1st1 <- stat_sig_s1st(s1st
#                      ,method="gtr"
#                      ,greater=0.01
#                      ,sigCri='either')

# S2: using the confidence intervals
s2_sig1 <- stat_sig_s2(s2
                       ,s2_conf_low
                       ,s2_conf_high
                       ,method='congtr'
                       ,greater=sig.cutoff
                       )

# S2: using greater than a given value
#s2_sig1 <- stat_sig_s2(s2
#                       ,s2_conf
#                       ,greater=0.02
#                       ,method='gtr')


##==============================================================================
## Barplots
##

# - along x axis, have the parameter names
# - along y axis, have first-order indices and total sensitivity indices
#   for the statistically significant parameters only
# - if parameter is not statistically significant, just leave out

par_indices <- seq(1,n_params)
ind_s1_sig <- which(s1st1$s1_sig==1)
ind_st_sig <- which(s1st1$st_sig==1)

icol_st <- 3
icol_s1 <- 11

pdf(paste(plotdir,'sobol_indicies.pdf',sep=''), width=8,height=5,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(1.4,.7,.15,.2))
ii <- 1
plot(c(ii,ii), c(0,s.out$T$original[ii]), col='azure4', lwd=1.5, type='l',
     xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, 57), ylim=c(0, 0.74))
for (ii in ind_st_sig) {lines(c(ii,ii), c(0,s.out$T$original[ii]), col='azure4', lwd=1.5)}
points(par_indices[ind_st_sig], s.out$T$original[ind_st_sig],
       col=rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), pch=16, cex=1.2)
points(par_indices[ind_s1_sig], s.out$S$original[ind_s1_sig],
       col=rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3]), pch=16, cex=1.2)
axis(1, at=par_indices, labels=parnames.sobol, cex.axis=.8, las=2)
axis(2, at=c(0,.1,.2,.3,.4,.5,.6,.7), labels=c('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'), cex.axis=.8, las=1)
mtext('Parameter', side=1, line=5.4, cex=.8)
mtext('Sensitivity index', side=2, line=2.4, cex=.8)
legend(1, 0.7, c('Total sensitivity', 'First-order'), pch=c(16,16), cex=.8, bty='n',
       col=c(rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3])))
dev.off()






##==============================================================================
## Setting up for radial plots
##

# Define groups for the variables and the color schemes
# Defining lists of the variables for each group

name_list1 <- list('All' = parnames.sobol)

# add Parameter symbols to plot
name_symbols <- c('ACT', expression('ACT'['carb']), 'VNV', 'NV', expression('e'^('NV')),
                  'LIFE', 'GYM', 'FERT', expression('e'^('fnBb')),
                  'dT2X', 'GLAC', 'J', 'n', 'Ws', expression('e'^('fD')), expression('Fwpa'['0']),
                  expression('Fwsa'['0']), expression('Fwga'['0']), expression('Fwca'['0']),
                  expression('Fmg'['0']), expression('Fmc'['0']), expression('Fmp'['0']),
                  expression('Fms'['0']), expression('Fwsi'['0']), expression('Xvolc'['0']),
                  expression('CAPd13C'['0']), expression('CAPd34S'['0']), expression('oxy'['570']),
                  expression('Gy'['570']), expression('Cy'['570']), expression('Ca'['570']),
                  expression('Ssy'['570']), expression('Spy'['570']), expression('dlsy'['570']),
                  expression('dlcy'['570']), expression('dlpy'['570']), expression('dlpa'['570']),
                  expression('dlgy'['570']), expression('dlga'['570']), expression('Rcy'['570']),
                  expression('Rca'['570']), expression('Rv'['570']), expression('Rg'['570']),
                  'Fob', 'COC', 'Ga', 'Ssa', 'Spa', 'ST', 'dlst', 'CT', 'dlct',
                  'kwpy', 'kwsy', 'kwgy', 'kwcy'
)


# defining list of colors for each group
col_list1 <- list("All"     = rgb(mycol[11,1],mycol[11,2],mycol[11,3]))

# using function to assign variables and colors based on group
s1st1 <- gp_name_col(name_list1
                     ,col_list1
                     ,s1st1)

s1st1$symbols <- name_symbols


plot.filename <- paste(plotdir,'sobol_spider',sep='')

plotRadCon(df=s1st1
           ,s2=s2
           ,plotS2=FALSE
           ,scaling = .36
           ,s2_sig=s2_sig1
           ,filename = plot.filename
           ,plotType = 'EPS'
           ,gpNameMult=15
           ,RingThick=0.1
           ,legLoc = "bottomcenter"
           ,cex = .6
           ,s1_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,st_col = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
           ,line_col = rgb(mycol[10,1],mycol[10,2],mycol[10,3])
           ,STthick = 0.5
           ,legFirLabs=c(.05,.77), legTotLabs=c(.05,.83), legSecLabs=c(.02,.05)
)

##
## Further analysis for the text:
##

# what are the highest first-order indices?
s1.sort <- s1st[rev(order(s1st[,'S1'])),1:4]
itmp <- which(s1.sort[,'S1'] > sig.cutoff & s1.sort[,'S1_conf_low']*s1.sort[,'S1_conf_high'] > 0)
s1.sort <- s1.sort[itmp,]
print('********************************')
print('significant first-order indices:')
print(s1.sort)
print('********************************')

# what are the highest total-order indices?
st.sort <- s1st[rev(order(s1st[,'ST'])),c(1,5:7)]
itmp <- which(st.sort[,'ST'] > sig.cutoff & st.sort[,'ST_conf_low']*st.sort[,'ST_conf_high'] > 0)
st.sort <- st.sort[itmp,]
print('********************************')
print('significant total-order indices:')
print(st.sort)
print('********************************')

# what are the highest second-order interaction indices?
s2.sort <- s2_table[rev(order(s2_table[,3])),]
itmp <- which(s2.sort[,'S2'] > sig.cutoff & s2.sort[,'S2_conf_low']*s2.sort[,'S2_conf_high'] > 0)
s2.sort <- s2.sort[itmp,]
print('********************************')
print('significant second-order indices:')
print(s2.sort)
print('********************************')

##==============================================================================



##==============================================================================
## End
##==============================================================================





##==============================================================================
## Try Morris?

alpha <- 0.1
perc_lower <- rep(0.5*alpha  , n_parameters)
perc_upper <- rep(1-0.5*alpha, n_parameters)

repetitions <- c(20, 40, 80)
levels <- c(20, 40, 80)

n_rep <- length(repetitions)
n_lev <- length(levels)

names_rep <- paste('r',repetitions[1],sep='')
if(n_rep>1) {for (rr in 2:n_rep) {names_rep <- c(names_rep, paste('r',repetitions[rr], sep=''))}}
names_lev <- paste('l',levels[1],sep='')
if(n_lev>1) {for (ll in 2:n_lev) {names_lev <- c(names_lev, paste('l',repetitions[ll], sep=''))}}

mom <- vector('list', n_rep); names(mom) <- names_rep
for (rr in names_rep) {
  mom[[rr]] <- vector('list', n_lev); names(mom[[rr]]) <- names_lev
}

## initialize other output
mu <- mu_star <- med <- med_star <- sigma <- stderr_mu <- mom_table <- names_signif <- mom

for (rr in 1:n_rep) {
  for (ll in 1:n_lev) {
    tbeg <- proc.time()
    grid_jump <- round(0.5*levels[ll])
    mom[[rr]][[ll]] <- morris(model=sensitivity_co2, factors=parnames_calib, r=repetitions[rr],
                              binf=perc_lower, bsup=perc_upper, scale=FALSE,
                              design=list(type="oat", levels=levels[ll], grid.jump=grid_jump),
                              par_fixed=par_fixed0, parnames_calib=parnames_calib,
                              parnames_fixed=parnames_fixed, age=age, ageN=ageN,
                              ind_const_calib=ind_const_calib, ind_time_calib=ind_time_calib,
                              ind_const_fixed=ind_const_fixed, ind_time_fixed=ind_time_fixed,
                              input=input,
                              ind_expected_time=ind_expected_time, ind_expected_const=ind_expected_const,
                              iteration_threshold=iteration_threshold, l_scaled=FALSE)
    tend <- proc.time()
    print(paste(length(mom[[rr]][[ll]]$y),' simulations took ',(tend-tbeg)[3]/60,' minutes',sep=''))
    mu[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, mean)
    mu_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) mean(abs(x)))
    med[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, median)
    med_star[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, function(x) median(abs(x)))
    sigma[[rr]][[ll]] <- apply(mom[[rr]][[ll]]$ee, 2, sd)
    stderr_mu[[rr]][[ll]] <- sigma[[rr]][[ll]]/sqrt(repetitions[rr])
    mom_table[[rr]][[ll]] <- cbind(mu[[rr]][[ll]], mu_star[[rr]][[ll]], sigma[[rr]][[ll]], stderr_mu[[rr]][[ll]])
  }
}
## Get the set of significant (mu_star > 2*stderr_mu) parameters for each of the
## simulations; take the intersection (union?) of the "okay" ones as the
## parameter set for calibration?
names_all <- NULL
for (rr in 1:n_rep) {
  for (ll in 1:n_lev) {
    names_signif[[rr]][[ll]] <- parnames_calib[which(mu_star[[rr]][[ll]] > 2*stderr_mu[[rr]][[ll]])]
    names_all <- union(names_all, names_signif[[rr]][[ll]])
  }
}

##==============================================================================



##==============================================================================


##==============================================================================
## Save progress
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
save.image(file=paste('geocarb_precalibration_',today,'.RData', sep=''))


##==============================================================================



##==============================================================================
## End
##==============================================================================
