##==============================================================================
## process_results.R
##
## read results, compute convergence diagnostics, thin chains
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================

setwd('~/codes/GEOCARB/R')
library('coda')
library('Hmisc')
source("compute_maxlag.R")
source("compute_grdiag.R")
source("pick_chains.R")
chains <- NULL



## [1] #########################################################################
##        dPpPUsOlUsn
appen <- "dPpPUsOlUsn"
datestamp <- "10Sep2019"
source("process_supp_chains.R");


## [2] #########################################################################
##        dPpPUsLlUsn
appen <- "dPpPUsLlUsn"
datestamp <- "11Sep2019"
source("process_supp_chains.R");


## [3] #########################################################################
##        dPpPUsRlUsn
appen <- "dPpPUsRlUsn"
datestamp <- "10Sep2019"
source("process_supp_chains.R");


## [4] #########################################################################
##        dPpPUsRlMsn
appen <- "dPpPUsRlMsn"
datestamp <- "10Sep2019"
source("process_supp_chains.R");


## [5] #########################################################################
##        dFpPUsRlUsn
appen <- "dFpPUsRlUsn"
datestamp <- "15Sep2019"
source("process_supp_chains.R");


## [6] #########################################################################
##        dFpPUsRlMsn
appen <- "dFpPUsRlMsn"
datestamp <- "10Sep2019"
source("process_supp_chains.R");


## [7] #########################################################################
##        dPpAUsRlMsn
## Paper results: use chains 1 and 3
appen <- "dPpAUsRlMsn"
datestamp <- "13Sep2019"
{source("process_supp_chains.R")};


## [8] #########################################################################
##        dFpAUsRlUsn
## Paper results: use chains 1, 3, 4 and 5
appen <- "dFpAUsRlUsn"
datestamp <- "18Sep2019"
{source("process_supp_chains.R")};


## [9] #########################################################################
##        dFpAUsRlMsn
appen <- "dFpAUsRlMsn"
datestamp <- "19Sep2019"
{source("process_supp_chains.R")};


## [10] #########################################################################
##        dFpAUsRlMnm
## Paper results: use chains 1, 2, 3 and 5
appen <- "dFpAUsRlMnm"
datestamp <- "18Sep2019"
{source("process_supp_chains.R")};


################################################################################

# saving for now to work on something else
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_chains <- paste('../output/chains_analysis_',today,'.rds', sep="")
saveRDS(chains, file=filename_chains)

# start back up again
chains <- readRDS("chains_analysis_30Sep2019.rds")

##
## compute ESS quantiles
##

quantiles_i_want <- c(.025, 0.05, 0.17, 0.25, 0.5, 0.75, 0.83, 0.95, 0.975)

ess <- vector('list', 2)
names(ess) <- c("ess", "ess_glac")
for (ee in names(ess)) {
  ess[[ee]] <- mat.or.vec(length(chains), length(quantiles_i_want))
  rownames(ess[[ee]]) <- names(chains)
  colnames(ess[[ee]]) <- quantiles_i_want
  for (experiment in names(chains)) {
    if (ncol(chains[[experiment]])==7) {
      idx_ess <- 5
    } else {
      idx_ess <- 10
    }
    if (ee=="ess") {
      ess[[ee]][experiment,] <- quantile(chains[[experiment]][,idx_ess], quantiles_i_want)
    } else {
      ess[[ee]][experiment,] <- quantile(chains[[experiment]][,idx_ess]*chains[[experiment]][,idx_ess+1], quantiles_i_want)
    }
  }
}

## Get the quantiles from Park and Royer 2011:
pr2011_dat <- read.csv('../input_data/ParkRoyer2011_Fig3_85varred.csv')
pr2011_cdf <- approxfun(pr2011_dat[,4], pr2011_dat[,1])
pr2011_pdf <- approxfun(pr2011_dat[,1], pr2011_dat[,3])
x_pr2011 <- pr2011_cdf(quantiles_i_want)
ess$ess <- rbind(x_pr2011, ess$ess)

## Plots

col_x <- c(8, 10.4, 13.3, 15.8, 18.5)
col_names <- c("Data", "Parameters", "fSR", "Likelihood", "Uncertainties")

table_experiments <- mat.or.vec(nr=nrow(ess$ess), length(col_names))
colnames(table_experiments) <- col_names
for (row in 1:nrow(ess$ess)) {
    chainname <- rownames(ess$ess)[row]
    if (chainname == "x_pr2011") {table_experiments[row,] <- c("PR2011", "PR2011", "PR2011", "unimodal", "log-normal")
    } else {
        if (substr(chainname, 2, 2)=="P") {table_experiments[row, "Data"] <- "PR2011"} else {table_experiments[row, "Data"] <- "F2017"}
        if (substr(chainname, 4, 5)=="PU") {table_experiments[row, "Parameters"] <- "PR2011"} else {table_experiments[row, "Parameters"] <- "all"}
        if (substr(chainname, 7, 7)=="R") {table_experiments[row, "fSR"] <- "DT2019"} else if (substr(chainname, 7, 7)=="L") {table_experiments[row, "fSR"] <- "L2018"} else {table_experiments[row, "fSR"] <- "PR2011"}
        if (substr(chainname, 9, 9)=="U") {table_experiments[row, "Likelihood"] <- "unimodal"} else {table_experiments[row, "Likelihood"] <- "mixture"}
        if (substr(chainname, 10, 11)=="sn") {table_experiments[row, "Uncertainties"] <- "skew-normal"} else {table_experiments[row, "Uncertainties"] <- "normal"}
    }
}

pdf(file='../figures/boxplot_ess.pdf', width=8, height=3.7, colormodel="cmyk", pointsize=11)
offset <- 0.045
yhgt <- offset*.75
experiment <- rownames(ess$ess)[1]
par(mfrow=c(1,1), mai=c(.7,.2,.2,.2))
plot(ess$ess[experiment,"0.5"], yhgt, xlim=c(0,22), ylim=c(0, .58), pch=16,
     xaxs='i', yaxs='i', yaxt='n', ylab='', xlab='', axes=FALSE)
grid()
axis(1, at=seq(0,9,2)); axis(1, at=seq(0,8,1), labels=rep("", 9), lwd=0.25)
mtext(expression(Delta*"T(2x) ["*degree*"C]"), side=1, line=2.2, adj=0.18)
for (experiment in rownames(ess$ess)) {
    row <- which(rownames(ess$ess)==experiment)
    points(ess$ess[experiment,"0.5"], yhgt, pch=16)
    arrows(x0=ess$ess[experiment,"0.05"], x1=ess$ess[experiment,"0.95"], y0=yhgt, y1=yhgt, angle=90, length=0.05, code=3)
    #text(8, yhgt, experiment, pos=4)
    #text(8, yhgt, table_experiments[row,], pos=4)
    for (i in 1:length(col_x)) {text(col_x[i], yhgt, table_experiments[row, i], pos=4)}
    yhgt <- yhgt + offset
}
yhgt <- yhgt + offset*0.25
for (i in 1:length(col_x)) {text(col_x[i], yhgt, colnames(table_experiments)[i], pos=4)}
lines(c(8, 22), c(yhgt-0.4*offset, yhgt-0.4*offset))
dev.off()



##
## compute model quantiles
##

library(gdata)
dat <- read.xls("../input_data/GEOCARB_output--both error envelopes_TW.xlsx", sheet="GEOCARB_output_INNER_BANDS")
model_quantiles_pr2011 <- cbind(dat[,"X0.025"], dat[,"CO2"], dat[,"X0.975"])
colnames(model_quantiles_pr2011) <- c('0.025','0.5','0.975')

model_experiment_quantiles <- vector('list', length(chains))
names(model_experiment_quantiles) <- names(chains)
for (ee in names(model_experiment_quantiles)) {
  if(substr(ee, 2, 2)=="P") {data_choice <- "PR2011"
  } else {data_choice <- "F2017"}
  if(substr(ee, 4, 4)=="P") {param_choice <- "PR2011_stdev"
  } else {param_choice <- "all_stdev"}
  if(substr(ee, 7, 7)=="O") {fSR_choice <- "PR2011"
  } else if(substr(ee, 7, 7)=="L") {fSR_choice <- "LENTON"
  } else if(substr(ee, 7, 7)=="R") {fSR_choice <- "ROYER"}
  if(substr(ee, 9, 9)=="U") {lhood_choice <- "unimodal"
  } else {lhood_choice <- "mixture"}
  dist <- substr(ee, 10, 11)
  source('model_setup.R')
  if(data_choice=="PR2011") {data_calib_pr2011 <- data_calib
  } else if(data_choice=="F2017") {data_calib_f2017 <- data_calib}
  model_out <- sapply(X=1:nrow(chains[[ee]]),
                FUN=function(k){model_forMCMC(par_calib=chains[[ee]][k,],
                                              par_fixed=par_fixed0,
                                              parnames_calib=parnames_calib,
                                              parnames_fixed=parnames_fixed,
                                              parnames_time=parnames_time,
                                              age=age,
                                              ageN=ageN,
                                              ind_const_calib=ind_const_calib,
                                              ind_time_calib=ind_time_calib,
                                              ind_const_fixed=ind_const_fixed,
                                              ind_time_fixed=ind_time_fixed,
                                              ind_expected_time=ind_expected_time,
                                              ind_expected_const=ind_expected_const,
                                              iteration_threshold=iteration_threshold,
                                              do_sample_tvq=DO_SAMPLE_TVQ,
                                              par_time_center=par_time_center,
                                              par_time_stdev=par_time_stdev)[,'co2']})
  model_experiment_quantiles[[ee]] <- mat.or.vec(nr=n_time, nc=3)
  for (t in 1:n_time) {
    model_experiment_quantiles[[ee]][t,] <- quantile(model_out[t,], c(.025, .5, .975))
  }
  colnames(model_experiment_quantiles[[ee]]) <- c('0.025', '0.5', '0.975')
}
model_experiment_quantiles$pr2011 <- model_quantiles_pr2011

# TODO - 5 row x 2 col figure of each experiment relative to the original PR2011
# Include the data set for each

pdf('../figures/model_experiments_vs_pr2011.pdf',width=7,height=9,colormodel='cmyk')
par(mfrow=c(5,2))
par(cex = 0.75)
par(mar = c(0, 0, 1.5, 0), oma = c(3.5, 5, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
panel <- 1
for (ee in names(chains)) {
    chainname <- ee
    if (chainname == "x_pr2011") {table_experiments[row,] <- c("PR2011", "PR2011", "PR2011", "unimodal", "log-normal")
    } else {
        if (substr(chainname, 2, 2)=="P") {data_choice <- "PR2011"} else {data_choice <- "F2017"}
        if (substr(chainname, 4, 5)=="PU") {param_choice <- "PR2011"} else {param_choice <- "all"}
        if (substr(chainname, 7, 7)=="R") {fSR_choice <- "DT2019"} else if (substr(chainname, 7, 7)=="L") {fSR_choice <- "L2018"} else {fSR_choice <- "PR2011"}
        if (substr(chainname, 9, 9)=="U") {likelihood_choice <- "unimodal"} else {likelihood_choice <- "mixture"}
        if (substr(chainname, 10, 11)=="sn") {kernel_choice <- "skew-normal"} else {kernel_choice <- "normal"}
    }
    title_tag <- paste(data_choice,"data,",param_choice,"parameters,",fSR_choice,"fSR,",likelihood_choice,"likelihood,",kernel_choice,"uncertainties")

    if (substr(ee, 2,2)=="P") {data_calib <- data_calib_pr2011
    } else {data_calib <- data_calib_f2017}
    plot(-time, log10(model_experiment_quantiles[[ee]][,'0.5']), type='l', xlim=c(-450,0), ylim=c(2,log10(6000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
    points(-data_calib$age, log10(data_calib$co2), pch=4, cex=0.4, lwd=.4)
    for (ii in 1:nrow(data_calib)) {
        arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
        if(substr(ee, 2,2)=="F") {arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)}
    }
    polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles$pr2011[,'0.025'],rev(model_experiment_quantiles$pr2011[,'0.975']))), col=rgb(.7,.2,.4,.5), border=NA)
    lines(-time, log10(model_experiment_quantiles$pr2011[,'0.5']), lwd=1.5, lty=5, col=rgb(.6,.2,.6))
    polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles[[ee]][,'0.025'],rev(model_experiment_quantiles[[ee]][,'0.975']))), col=rgb(.2,.6,.6,.5), border=NA)
    lines(-time, log10(model_experiment_quantiles[[ee]][,'0.5']), lwd=1.5, col=rgb(.2,.6,.75))
    #mtext(ee, side=3, cex=0.8)
    mtext(title_tag, side=3, cex=0.8)
    if (panel %in% c(9,10)) {
        axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
        mtext('Time [Myr ago]', side=1, line=2.1, cex=0.8)
    }
    ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))
    if (panel %in% c(1,3,5,7,9)) {
        axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
        axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
        mtext(expression('CO'[2]*' [ppmv]'), side=2, line=3, cex=0.8)
    }
    if (panel==1) {
        legend(-390, log10(6600), c('95% range (Royer [2014])', '95% range (this work)'), pch=c(15,15), col=c(rgb(.7,.2,.4,.5), rgb(.2,.6,.6,.5)), cex=.8, bty='n')
        legend(-180, log10(6600), c('',''), pch=c(NA,4), col=c('black','black'), cex=.8, bty='n')
        legend(-180, log10(6600), c('Median','Data'), pch=c('-',NA), col=c('black','black'), cex=.8, bty='n')
    }
    minor.tick(nx=5, ny=0, tick.ratio=0.5)
    panel <- panel+1
}
dev.off()

##==============================================================================



if(FALSE) {
##
## try to have a menu to pick which two simulation sets to plot
##

library(manipulate)

par(mfrow=c(2,1), mai=c(.8,.8,.1,.1))
manipulate(
    {
        runname <- ""
        if (dataset=="dF") {data_calib <- data_calib_f2017} else if (dataset=="dP") {data_calib <- data_calib_pr2011}
        runname <- paste(runname, dataset, sep="")
        runname <- paste(runname, parameters, sep="")
        runname <- paste(runname, seafloor_spreading, sep="")
        runname <- paste(runname, likelihood, sep="")
        runname <- paste(runname, kernels, sep="")
        exp1 <- "pr2011"
        exp2 <- runname; print(runname)
        if ((substr(exp1, 2,2)=="F") | (substr(exp2, 2,2)=="F")) {data_calib <- data_calib_f2017} else {data_calib <- data_calib_pr2011}
        plot(-time, log10(model_experiment_quantiles[[exp1]][,'0.5']), type='l', xlim=c(-450,0), ylim=c(2,log10(6000)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n')
        points(-data_calib$age, log10(data_calib$co2), pch=4, cex=0.4, lwd=.4)
        for (ii in 1:nrow(data_calib)) {
            arrows(-data_calib$age[ii], log10(data_calib$co2_low[ii]), -data_calib$age[ii], log10(data_calib$co2_high[ii]), length=0.02, angle=90, code=3, lwd=0.5)
            if(nrow(data_calib)==nrow(data_calib_f2017)) {arrows(-data_calib$age_old[ii], log10(data_calib$co2[ii]), -data_calib$age_young[ii], log10(data_calib$co2[ii]), length=0.02, angle=90, code=3, lwd=0.5)}
        }
        polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles[[exp2]][,'0.025'],rev(model_experiment_quantiles[[exp2]][,'0.975']))), col=rgb(.7,.2,.4,.5), border=NA)
        lines(-time, log10(model_experiment_quantiles[[exp2]][,'0.5']), lwd=1.5, lty=5, col=rgb(.6,.2,.6))
        polygon(-c(time,rev(time)), log10(c(model_experiment_quantiles[[exp1]][,'0.025'],rev(model_experiment_quantiles[[exp1]][,'0.975']))), col=rgb(.2,.6,.6,.5), border=NA)
        lines(-time, log10(model_experiment_quantiles[[exp1]][,'0.5']), lwd=1.5, col=rgb(.2,.6,.75))
        mtext(ee, side=3, cex=0.8)
        axis(1, at=seq(-400,0,100), labels=c('400','300','200','100','0'), cex.axis=1.1)
        mtext('Time [Myr ago]', side=1, line=2.1, cex=0.8)
        ticks=log10(c(seq(10,100,10),seq(200,1000,100),seq(2000,10000,1000)))

        axis(2, at=ticks, labels=rep('',length(ticks)), cex.axis=1.1)
        axis(2, at=log10(c(10,30,100,300,1000,3000)), labels=c('10','30','100','300','1000','3000'), cex.axis=1.1, las=1)
        mtext(expression('CO'[2]*' [ppmv]'), side=2, line=3, cex=0.8)


        legend(-390, log10(6600), c(paste('95% range',exp2), paste('95% range',exp1)), pch=c(15,15), col=c(rgb(.7,.2,.4,.5), rgb(.2,.6,.6,.5)), cex=.8, bty='n')
        legend(-180, log10(6600), c('',''), pch=c(NA,4), col=c('black','black'), cex=.8, bty='n')
        legend(-180, log10(6600), c('Median','Data'), pch=c('-',NA), col=c('black','black'), cex=.8, bty='n')

        minor.tick(nx=5, ny=0, tick.ratio=0.5)},
    dataset = picker("Park and Royer [2011]" = "dP", "Foster et al [2017]" = "dF"),
    parameters = picker("Park and Royer [2011]" = "pPU", "All" = "pAU"),
    seafloor_spreading = picker("Original" = "sO", "Lenton et al [2018]" = "sL", "Domeier and Torsvik [2019]" = "sR"),
    likelihood = picker("Unimodal" = "lU", "Mixture" = "lM"),
    kernels = picker("Skew-normal" = "sn", "Normal" = "nm"))

    #exp1 = picker("dPpPUsOlUsn", "dPpPUsLlUsn", "dPpPUsRlUsn", "dPpPUsRlMsn", "dFpPUsRlMsn", "dFpPUsRlUsn", "dPpAUsRlMsn", "dFpAUsRlMsn", "dFpAUsRlMnm", "pr2011"),
    #exp2 = picker("dPpPUsOlUsn", "dPpPUsLlUsn", "dPpPUsRlUsn", "dPpPUsRlMsn", "dFpPUsRlMsn", "dFpPUsRlUsn", "dPpAUsRlMsn", "dFpAUsRlMsn", "dFpAUsRlMnm", "pr2011"))

}


##==============================================================================
## END
##==============================================================================
