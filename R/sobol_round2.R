icol_st##==============================================================================
## sensitivity_geocarb.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())


n_sample <- 1e5
alpha <- 0.34


setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')



##==============================================================================
## Previous Sobol results
##=======================

# Read previous results
s.out <- readRDS('../output/sobol_alpha34_18Jan2018.rds')

# Whicih parameters to have a closer look at?
ind.sample <- which(s.out$T[,1] > 0.05)
ind.fix <- which(s.out$T[,1] <= 0.05)

# Get the parameter names
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')


# Create a parameters file for this calibration/sensitivity analysis
calib <- read.csv('../input_data/GEOCARB_input_summaries_calib_all.csv')

for (k in 1:length(ind.fix)) {
  row = which(calib$parameter==parnames_calib[ind.fix[k]])
  calib$calib[row] <- 0
}

# Write the revised file
write.csv(x=calib, file='../input_data/GEOCARB_input_summaries_calib_newSobol.csv', row.names=FALSE)

##==============================================================================


##==============================================================================
## Model parameters and setup
##===========================

# Read parameter information, set up the calibration parameters
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib_newSobol.csv'
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
alpha <- 1-alpha
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

tbeg <- proc.time()

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

tend <- proc.time()
print(paste('LHS precalibration simulations and filtering took ',tend[3]-tbeg[3],' seconds', sep=''))

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
filename_out <- paste('../output/geocarb_precalibration_parameters_alpha',100*alpha,'_',today,'.csv', sep='')
write.csv(rbind(parameters_good, bandwidths), file=filename_out, row.names=FALSE)
##==============================================================================


##==============================================================================
## Function for sampling from KDEs (inflating a LHS?)
##===================================================

## Read KDE results file, separate into parameters and the bandwidths
filename_in <- '../output/geocarb_precalibration_parameters_alpha34_04Feb2018.csv'
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

print(paste('Sobol simulations took ',t.out[3],' seconds', sep=''))

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.sobol <- paste('../output/sobol_alpha',100*alpha,'_',today,'.rds', sep='')
saveRDS(s.out, filename.sobol)

## Check convergence - is maximum confidence interval width < 10% of the highest
## total order index?
max_sens_ind <- 0.1*max(s.out$T$original)
max_conf_int <- max(max(s.out$S$`max. c.i.` - s.out$S$`min. c.i.`), max(s.out$T$`max. c.i.` - s.out$T$`min. c.i.`))
print(paste('0.1 * max. sensitivity index=',max_sens_ind,' // max. conf int=',max_conf_int,sep=''))
##==============================================================================



##==============================================================================
## Write indices, results we'd need to file

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output/geocarb_sobol-1-tot_alpha',100*alpha,'_',today,'.txt',sep='')
file.sobolout2 <- paste('../output/geocarb_sobol-2_alpha',100*alpha,'_',today,'.txt',sep='')

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
n_params <- 27

# Set Sobol indices file names
Sobol_file_1 <- "../output/geocarb_sobol-1-tot_alpha34_04Feb2018.txt"
Sobol_file_2 <- "../output/geocarb_sobol-2_alpha34_04Feb2018.txt"
s.out <- readRDS('../output/sobol_alpha34_04Feb2018.rds')

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

pdf(paste(plotdir,'sobol_indices_new.pdf',sep=''), width=8,height=5,colormodel='cmyk')
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

# what is the total variance in output accounted for by first- and second-order
# interations?
print('********************************')
print(sum(s1st1$S1[which(s1st1$s1_sig==1)]) + sum(s2 * s2_sig1, na.rm=TRUE))
print('********************************')
##==============================================================================





##==============================================================================
## Write calibration input files
##==============================


##==========
# Full set of total sensitivity index > 0.05 from original Sobol'
rv <- file.copy(from=filename.calibinput, to='../input_data/GEOCARB_input_summaries_calib_ST05.csv', overwrite=TRUE)
##==========


##==========
# Only the parameters with significant total index after second Sobol'
ind.fix <- which(s1st1$st_sig==0)


#TODO
##==========


##==========
# Only the parameters with significant first- or second-order index after second Sobol'
ind_s1_sig <- which(s1st1$s1_sig==1)

ind_s2_sig <- NULL
for (p in 1:ncol(s2_sig1)) {
    ind_s2_sig <- c(ind_s2_sig, which(s2_sig1[,p]==1))
}
ind_s2_sig <- unique(ind_s2_sig)

ind_sig <- unique(c(ind_s1_sig, ind_s2_sig))

#TODO


ind.fix <- which(s.out$T[,1] <= 0.05)

# Get the parameter names
filename.calibinput <- '../input_data/GEOCARB_input_summaries_calib.csv'
source('GEOCARB-2014_parameterSetup.R')

# Create a parameters file for this calibration/sensitivity analysis
calib <- read.csv('../input_data/GEOCARB_input_summaries_calib_all.csv')

calib$calib <- 0

for (k in 1:length(ind_sig)) {
  row = which(calib$parameter==parnames_calib[ind_sig[k]])
  calib$calib[row] <- 1
}

# Write the revised file
write.csv(x=calib, file='../input_data/GEOCARB_input_summaries_calib_S1S2.csv', row.names=FALSE)
##==========

##==============================================================================
## End
##==============================================================================
