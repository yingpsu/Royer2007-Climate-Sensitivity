##==============================================================================
## plotting_sobol1T.R
##
## Plotting of...
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================


setwd('~/codes/GEOCARB/R')
Sobol_file_1 <- "../output/geocarb_sobol-1-tot_alpha0_sensNS_NS-n40K-bs10K_25Jun2018.txt"

n_params <- 56
plotdir <- '../figures/'

##==============================================================================

## some analysis

## Libraries
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

pdf(paste(plotdir,'sobol_indices.pdf',sep=''), width=8,height=5,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(1.4,.7,.15,.2))
ii <- 1
plot(c(ii,ii), c(0, s1st$ST[ii]), col='azure4', lwd=1.5, type='l',
     xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, 57), ylim=c(0, 0.74))
for (ii in ind_st_sig) {lines(c(ii,ii), c(0,s1st$ST[ii]), col='azure4', lwd=1.5)}
points(par_indices[ind_st_sig], s1st$ST[ind_st_sig],
       col=rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), pch=16, cex=1.2)
points(par_indices[ind_s1_sig], s1st$S1[ind_s1_sig],
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
## Plot and analyze the whole bunch

setwd('/Users/tony/codes/GEOCARB/R')
plotdir <- '../figures/'

## Libraries
library(RColorBrewer) # good color palettes
library(graphics)     # used when plotting polygons
library(plotrix)      # used when plotting circles

## Functions in other files
source('sobol_functions.R')
source('colorblindPalette.R')

# Set number of parameters being analyzed
n_params <- 56

# experiments
sobol_exp <- c('a50','a34','a10','a05'); n_experiment <- length(sobol_exp)
sobol_files <- vector('list', n_experiment); names(sobol_files) <- sobol_exp
s1st <- vector('list', n_experiment); names(s1st) <- sobol_exp
s1st1 <- vector('list', n_experiment); names(s1st1) <- sobol_exp
ind_s1_sig <- vector('list', n_experiment); names(ind_s1_sig) <- sobol_exp
ind_st_sig <- vector('list', n_experiment); names(ind_st_sig) <- sobol_exp

## load relevant output files

# Set Sobol indices file names

sobol_files$a50 <- "../output/geocarb_sobol-1-tot_alpha50_17Jan2018.txt"
sobol_files$a34 <- "../output/geocarb_sobol-1-tot_alpha34_18Jan2018.txt"
sobol_files$a10 <- "../output/geocarb_sobol-1-tot_alpha10_18Jan2018.txt"
sobol_files$a05 <- "../output/geocarb_sobol-1-tot_alpha05_18Jan2018.txt"

## Import data from sensitivity analysis
# First- and total-order indices
for (aa in sobol_exp) {
  s1st[[aa]] <- read.csv(sobol_files[[aa]],
                         sep=' ',
                         header=TRUE,
                         nrows = n_params,
                         as.is=c(TRUE,rep(FALSE,5)))
}
parnames.sobol <- s1st[[1]][,1]

## Determine which indices are statistically significant

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
for (aa in sobol_exp) {
  s1st1[[aa]] <- stat_sig_s1st(s1st[[aa]]
                               ,method="congtr"
                               ,greater=sig.cutoff
                               ,sigCri='either')
}

# S1 & ST: using greater than a given value
#for (aa in sobol_exp) {
#  s1st1[[aa]] <- stat_sig_s1st(s1st[[aa]]
#                               ,method="gtr"
#                               ,greater=sig.cutoff
#                               ,sigCri='either')
#}


for (aa in sobol_exp) {
  ind_s1_sig[[aa]] <- which(s1st1[[aa]]$s1_sig==1)
  ind_st_sig[[aa]] <- which(s1st1[[aa]]$st_sig==1)
}

st_sig <- mat.or.vec(n_params, n_experiment)
for (ii in 1:n_experiment) {
  st_sig[s1st1[[ii]]$ST_conf_low > sig.cutoff, ii] <- 1
}

cutoff_test <- seq(0.01, 0.15, by=0.01)
n_test <- length(cutoff_test)
st_sig_exp <- mat.or.vec(n_params, n_test)
cutoff_matches <- rep(FALSE, n_test)
for (cc in 1:n_test) {
  st_sig_exp[s1st1$a10$ST > cutoff_test[cc], cc] <- 1
  if(all(st_sig_exp[,cc]==s1st1$a50$st_sig)) {cutoff_matches[cc] <- TRUE}
}

wide.cutoff <- cutoff_test[which(cutoff_matches)[1]]
print(paste('cutoff we are using: ',wide.cutoff,sep=''))


# returns 0.06 and 0.07 for exact matches, so use cut-off of 0.05% of total
# variance as cut-off for "importance"

## Barplots

par_indices <- seq(1,n_params)

icol_st <- 3
icol_s1 <- 11


pdf(paste(plotdir,'sobol_indices_all.pdf',sep=''), width=8,height=10,colormodel='cmyk')
par(mfrow=c(3,1), mai=c(1.2,.7,.15,.2))
aa='a50'  # alpha=0.50
ii <- 1
plot(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5, type='l',
     xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, 57), ylim=c(0, 1.05))
for (ii in ind_st_sig[[aa]]) {lines(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5)}
points(par_indices[ind_st_sig[[aa]]], s1st1[[aa]]$ST[ind_st_sig[[aa]]],
       col=rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), pch=16, cex=1.2)
points(par_indices[ind_s1_sig[[aa]]], s1st1[[aa]]$S1[ind_s1_sig[[aa]]],
       col=rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3]), pch=16, cex=1.2)
axis(1, at=par_indices, labels=parnames.sobol, cex.axis=1, las=2)
axis(2, at=seq(0,1,.1), labels=c('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'), cex.axis=1, las=1)
mtext('Parameter', side=1, line=7, cex=.8)
mtext('Sensitivity index', side=2, line=3, cex=.8)
legend(1, 1, c('Total sensitivity', 'First-order'), pch=c(16,16), cex=1, bty='n',
       col=c(rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3])))

par(mai=c(1.2,.7,.15,.2))
aa='a10'  # alpha=0.10
ii <- 1
plot(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5, type='l',
     xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, 57), ylim=c(0, 1.05))
for (ii in ind_st_sig[[aa]]) {lines(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5)}
points(par_indices[ind_st_sig[[aa]]], s1st1[[aa]]$ST[ind_st_sig[[aa]]],
       col=rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), pch=16, cex=1.2)
points(par_indices[ind_s1_sig[[aa]]], s1st1[[aa]]$S1[ind_s1_sig[[aa]]],
       col=rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3]), pch=16, cex=1.2)
lines(c(0, n_params+1), c(wide.cutoff,wide.cutoff), lwd=1.5, type='l', lty=2, col='azure4')
axis(1, at=par_indices, labels=parnames.sobol, cex.axis=1, las=2)
axis(2, at=seq(0,1,.1), labels=c('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'), cex.axis=1, las=1)
mtext('Parameter', side=1, line=7, cex=.8)
mtext('Sensitivity index', side=2, line=3, cex=.8)

par(mai=c(1.2,.7,.15,.2))
aa='a05'  # alpha=0.05
ii <- 1
plot(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5, type='l',
     xaxt='n', yaxt='n', xlab='', ylab='', xaxs='i', yaxs='i', xlim=c(0, 57), ylim=c(0, 1.05))
for (ii in ind_st_sig[[aa]]) {lines(c(ii,ii), c(0,s1st1[[aa]]$ST[ii]), col='azure4', lwd=1.5)}
points(par_indices[ind_st_sig[[aa]]], s1st1[[aa]]$ST[ind_st_sig[[aa]]],
       col=rgb(mycol[icol_st,1],mycol[icol_st,2],mycol[icol_st,3]), pch=16, cex=1.2)
points(par_indices[ind_s1_sig[[aa]]], s1st1[[aa]]$S1[ind_s1_sig[[aa]]],
       col=rgb(mycol[icol_s1,1],mycol[icol_s1,2],mycol[icol_s1,3]), pch=16, cex=1.2)
lines(c(0, n_params+1), c(wide.cutoff,wide.cutoff), lwd=1.5, type='l', lty=2, col='azure4')
axis(1, at=par_indices, labels=parnames.sobol, cex.axis=1, las=2)
axis(2, at=seq(0,1,.1), labels=c('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'), cex.axis=1, las=1)
mtext('Parameter', side=1, line=7, cex=.8)
mtext('Sensitivity index', side=2, line=3, cex=.8)

dev.off()





##==============================================================================




##==============================================================================
## End
##==============================================================================
