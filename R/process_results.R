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
## Paper results: use chains 1, 2, 4, 5
appen <- "dPpAUsRlMsn"
datestamp <- "03Oct2019"
{source("process_supp_chains.R")};


## [8] #########################################################################
##        dFpAUsRlUsn
## Paper results: use chains 1, 3, 4, 5
appen <- "dFpAUsRlUsn"
datestamp <- "03Oct2019"
{source("process_supp_chains.R")};


## [9] #########################################################################
##        dFpAUsRlMsn
appen <- "dFpAUsRlMsn"
datestamp <- "04Oct2019"
{source("process_supp_chains.R")};


## [10] #########################################################################
##        dFpAUsRlMnm
appen <- "dFpAUsRlMnm"
datestamp <- "03OCt2019"
{source("process_supp_chains.R")};


################################################################################

# saving for now to work on something else
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename_chains <- paste('../output/chains_analysis_',today,'.rds', sep="")
saveRDS(chains, file=filename_chains)

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
