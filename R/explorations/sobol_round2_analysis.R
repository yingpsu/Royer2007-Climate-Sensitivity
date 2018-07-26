##==============================================================================
## sensitivity_geocarb.R
##
## Sensitivity experiment with GEOCARB model (Foster et al 2017 version)
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


## Clear workspace
rm(list=ls())

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')


##==============================================================================

## load relevant output files

setwd('/Users/tony/codes/Royer2007-Climate-Sensitivity/R')
plotdir <- '../figures/'

# Set number of parameters being analyzed
n_params <- 13

# Set Sobol indices file names
Sobol_file_1 <- "../output/geocarb_round2_sobol-1-tot_alpha34_05Feb2018.txt"
Sobol_file_2 <- "../output/geocarb_round2_sobol-2_alpha34_05Feb2018.txt"
s.out <- readRDS('../output/sobol_round2_alpha34_05Feb2018.rds')

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
