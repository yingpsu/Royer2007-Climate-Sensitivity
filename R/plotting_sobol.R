##==============================================================================
## plotting_sobol.R
##
## Plotting of...
##
## Questions?  Tony Wong (anthony.e.wong@colorado.edu)
##==============================================================================

#rm(list=ls())

#setwd('~/codes/GEOCARB/R')

filename.calibinput <- "../input_data/GEOCARB_input_summaries_calib_tvq_all.csv"

Sobol_file_1 <- "../output/geocarb_sobol-1-tot_sensNS_N47751-BS5k_10Oct2018.txt"
Sobol_file_2 <- "../output/geocarb_sobol-2_sensNS_N47751-BS5k_10Oct2018.txt"

#-rw-r--r--   1 tony  staff   7.8K Sep 28 08:38 geocarb_sobol-1-tot_sensNS_30kN-5kBS_28Sep2018.txt
#-rw-r--r--   1 tony  staff   148K Sep 28 08:38 geocarb_sobol-2_sensNS_30kN-5kBS_28Sep2018.txt
#-rw-r--r--   1 tony  staff    16M Sep 28 08:38 sobol_sensNS_30kN-5kBS_28Sep2018.rds

n_params <- 68
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

sig.cutoff <- 0.05

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
source('GEOCARB-2014_parameterSetup_tvq.R')

# yields parnames_calib, which are the sensitive parameters
ind_sensit <- which((parnames.sobol %in% parnames_calib) & (s1st1$sig>0))
ind_insens <- 1:n_params; ind_insens <- ind_insens[-ind_sensit]

name_list1 <- list('Sensitive' = parnames.sobol[ind_sensit]
               ,'Insensitive' = parnames.sobol[ind_insens]
               )

# add Parameter symbols to plot

name_symbols <- c('ACT', expression('ACT'['carb']), 'VNV', 'NV', expression('e'^'NV'),
                  'LIFE', 'GYM', 'FERT', expression('e'^'fnBb'),
                  expression(Delta*'T(2x)'), 'GLAC', 'J', 'n', 'Ws', expression('e'^'fD'), expression('Fwpa'['0']),
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
                  'kwpy', 'kwsy', 'kwgy', 'kwcy',
                  "Sr", "d13C", "d34S", "fR", "fL", "fA", "fD", "fAw_fA", "RT", "GEOG", "fSR", "fC"
)

new_name_symbols <- c(name_symbols[ind_sensit], name_symbols[ind_insens])

# defining list of colors for each group
col_list1 <- list("Sensitive"     = 'black'
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

s1st1_rearr <- rbind(s1st1[ind_sensit,], s1st1[ind_insens,])
s1st1_rearr$symbols <- new_name_symbols # it was right already!
s2_rearr <- rbind(s2[ind_sensit,],s2[ind_insens,]) # rearrange the rows...
s2_rearr <- cbind(s2_rearr[,ind_sensit],s2_rearr[,ind_insens]) # ... and the columns
s2_sig1_rearr <- rbind(s2_sig1[ind_sensit,], s2_sig1[ind_insens,])
s2_sig1_rearr <- cbind(s2_sig1_rearr[,ind_sensit], s2_sig1_rearr[,ind_insens])

##==============================================================================
## Plot only the sensitive parameters

ind_keep <- which(s1st1_rearr$sig==1)
s1st1_sens <- s1st1_rearr[ind_keep,]
s2_sens <- s2_rearr[ind_keep,ind_keep]
s2_sig1_sens <- s2_sig1_rearr[ind_keep,ind_keep]
s1st1_sens$symbols <- new_name_symbols[ind_keep]


plot.filename <- paste(plotdir,'sobol_spider_tvq_sensOnly',sep='')

plotRadCon(df=s1st1_sens
           ,s2=s2_sens
           ,plotS2=TRUE
           ,radSc = 2
           ,scaling=.5
           ,widthSc = 0.4
           ,s2_sig=s2_sig1_sens
           ,filename = plot.filename
           ,plotType = 'EPS'
           ,gpNameMult=100
           ,varNameMult=1.34
           ,RingThick=0.15
           ,legLoc = "bottomcenter"
           ,cex = 1
           ,rt_names = 0
           ,s1_col = 'coral'
           ,st_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,line_col = 'slategray2'
           ,STthick = 0.5
           ,legFirLabs=c(.05,.15), legTotLabs=c(.05,.75), legSecLabs=c(.02,.10)
)


##==============================================================================
## Some numbers to report

#We find that parameters related to the paleoclimate time series, in particular
#with respect to estimated weathering fluxes (RT, fAw_fA and fD) and seafloor
#spreading and degassing (fSR and fC), are sensitive (Fig. 3). The first- and
#second-order effects of uncertainty in all 12 time series required to force the
# GEOCARB model is XXX% of the total variance in model output. From this, we
#conclude that previous results using the various incarnations of the GEOCARB
#model are still indeed valid.


s2_time <- as.matrix(s2[,parnames_time])
# the s2 one will work because the time parameters are all at the end
variance_from_time <- sum(s1st1$S1[s1st1$Parameter %in% parnames_time]) +
                      sum(s2_time[which(s2_time>0)])
print(paste('Total first- and second-order variance contribution from time arrays (and their interactions with all other parameters) =',variance_from_time))


# Chemical weathering
parnames_to_add <- c('ACT','ACTcarb','VNV', 'NV', 'exp_NV')
s2_chem <- s2[parnames_to_add,]
variance_from_chemicalweathering <- 0
for (name in parnames_to_add) {
    variance_from_chemicalweathering <- variance_from_chemicalweathering + s1st1[match(name,s1st1$Parameter),'S1']
    # going by row works because the chem weathering params are all at the beginning
    ind_add <- which(s2_chem[name,]>0)
    if(length(ind_add)>0) {variance_from_chemicalweathering <- variance_from_chemicalweathering + sum(s2_chem[name,ind_add], na.rm=TRUE)}
}
print(paste('Total first- and second-order variance contribution from chemical weathering (and their interactions with all other parameters) =',variance_from_chemicalweathering))


# Plant-assisted weathering
parnames_to_add <- c('LIFE','GYM','FERT', 'exp_fnBb')
mindex <- min(match(parnames_to_add, parnames_calib))
s2_plant <- s2[parnames_to_add,]
variance_from_plantweathering <- 0
for (name in parnames_to_add) {
  variance_from_plantweathering <- variance_from_plantweathering + s1st1[match(name,s1st1$Parameter),'S1']
  # need to add the rest of this parameter's row
  ind_add <- which(s2_plant[name,]>0)
  if(length(ind_add)>0) {variance_from_plantweathering <- variance_from_plantweathering + sum(s2_plant[name,ind_add], na.rm=TRUE)}
  # and need to add this parameter's column **up through the first one in the group**
  ind_add <- which(s2_plant[1:mindex,name]>0)
  if(length(ind_add)>0) {variance_from_plantweathering <- variance_from_plantweathering + sum(s2_plant[ind_add,name], na.rm=TRUE)}
}
print(paste('Total first- and second-order variance contribution from plant-assisted weathering (and their interactions with all other parameters) =',variance_from_plantweathering))


# second-order sensitivity indices that are significant
print(paste('2nd order sensitivity index between deltaT2X and GYM =',s2_sens['GYM','deltaT2X']))
print(paste('2nd order sensitivity index between GYM and ACT =',s2_sens['ACT','GYM']))

##==============================================================================
## Radial convergence plot

plot.filename <- paste(plotdir,'sobol_spider_tvq',sep='')

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


##==============================================================================
## End
##==============================================================================
