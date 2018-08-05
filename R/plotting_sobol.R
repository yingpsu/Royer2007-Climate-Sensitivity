##==============================================================================
## plotting_sobol.R
##
## Plotting of...
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

rm(list=ls())

setwd('~/codes/GEOCARB/R')

Sobol_file_1 <- "../output/geocarb_sobol-1-tot_alpha0_sensNS_n30K-bs10K_03Aug2018.txt"
Sobol_file_2 <- "../output/geocarb_sobol-2_alpha0_sensNS_n30K-bs10K_03Aug2018.txt"
filename.calibinput <- "../input_data/GEOCARB_input_summaries_calib_sig18.csv"

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
colnames(s2) <- rownames(s2) <- parnames.sobol

# Convert confidence intervals to upper-triangular matrix
s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

s2_conf_low <- as.data.frame(s2_conf_low)
s2_conf_high <- as.data.frame(s2_conf_high)
colnames(s2_conf_low) <- rownames(s2_conf_low) <- parnames.sobol
colnames(s2_conf_high) <- rownames(s2_conf_high) <- parnames.sobol



##==============================================================================
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
## Get sensitive parameter names
source('GEOCARB-2014_parameterSetup.R')

# yields parnames_calib, which are the sensitive parameters
ind_sensit <- which(parnames.sobol %in% parnames_calib)
ind_insens <- 1:n_params; ind_insens <- ind_insens[-ind_sensit]

name_list1 <- list('Sensitive' = parnames.sobol[ind_sensit]
               ,'Insensitive' = parnames.sobol[ind_insens]
               )

# add Parameter symbols to plot


#
# TODO -- THESE NEED TO BE REARRANGED SO THAT THE SENSTIVIE PARAMETERS ARE
#         SHUFFLED DOWN TO THE FRONT, AND THE INSENSITIVE ONES TO THE BACK.
#
# NEEDS TO MATCH THE ORDER c(parnames.sobol[ind_sensit], parnames.sobol[ind_insens])
#

name_symbols <- c('ACT', expression('ACT'['carb']), 'VNV', 'NV', expression('e'^('NV')),
                  'LIFE', 'GYM', 'FERT', expression('e'^('fnBb')),
                  'dT2X', 'GLAC', 'J', 'n', 'Ws', expression('e'^('fD')), expression('Fwpa'['0']),
                  expression('Fwsa'['0']), expression('Fwga'['0']), expression('Fwca'['0']),
                  expression('Fmg'['0']), expression('Fmc'['0']), expression('Fmp'['0']),
                  expression('Fms'['0']), expression('Fwsi'['0']), expression('Xvolc'['0']),
                  expression('FRd13C'['0']), expression('FRd34S'['0']), expression('oxy'['570']),
                  expression('Gy'['570']), expression('Cy'['570']), expression('Ca'['570']),
                  expression('Ssy'['570']), expression('Spy'['570']), expression('dlsy'['570']),
                  expression('dlcy'['570']), expression('dlpy'['570']), expression('dlpa'['570']),
                  expression('dlgy'['570']), expression('dlga'['570']), expression('Rcy'['570']),
                  expression('Rca'['570']), expression('Rv'['570']), expression('Rg'['570']),
                  'Fob', 'COC', 'Ga', 'Ssa', 'Spa', 'ST', 'dlst', 'CT', 'dlct',
                  'kwpy', 'kwsy', 'kwgy', 'kwcy'
)

new_name_symbols <- c(name_symbols[ind_sensit], name_symbols[ind_insens])

# defining list of colors for each group
col_list1 <- list("Sensitive"     = rgb(mycol[11,1],mycol[11,2],mycol[11,3])
                  ,"Insensitive" = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
                  )

# using function to assign variables and colors based on group
s1st1 <- gp_name_col(name_list1
                     ,col_list1
                     ,s1st1)

s1st1$symbols <- new_name_symbols


# Things to change order of to group by sensitivity:
# 1) s1st1
# 2) s2
# 3) s2_sig1
# 4)

s1st1_rearr <- rbind(s1st1[ind_sensit,], s1st1[ind_insens,])
s1st1_rearr$symbols <- new_name_symbols # it was right already!
s2_rearr <- rbind(s2[ind_sensit,],s2[ind_insens,]) # rearrange the rows...
s2_rearr <- cbind(s2_rearr[,ind_sensit],s2_rearr[,ind_insens]) # ... and the columns
s2_sig1_rearr <- rbind(s2_sig1[ind_sensit,], s2_sig1[ind_insens,])
s2_sig1_rearr <- cbind(s2_sig1_rearr[,ind_sensit], s2_sig1_rearr[,ind_insens])

##==============================================================================
## Radial convergence plot (aka sobol spider diagram)
##

plot.filename <- paste(plotdir,'sobol_spider',sep='')

plotRadCon(df=s1st1_rearr
           ,s2=s2_rearr
           ,plotS2=TRUE
           ,radSc = 2
           ,scaling=.3
           ,widthSc = 0.5
           ,s2_sig=s2_sig1_rearr
           ,filename = plot.filename
           ,plotType = 'EPS'
           ,gpNameMult=1.7
           ,varNameMult=1.34
           ,RingThick=0.17
           ,legLoc = "bottomcenter"
           ,cex = .75
           ,s1_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,st_col = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
           ,line_col = rgb(mycol[10,1],mycol[10,2],mycol[10,3])
           ,STthick = 0.5
           ,legFirLabs=c(.05,.25), legTotLabs=c(.05,.75), legSecLabs=c(.02,.07)
)


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
