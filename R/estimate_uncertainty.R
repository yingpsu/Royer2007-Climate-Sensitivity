##==============================================================================
## estimate_uncertainty.R
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

n_time <- 58
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_25Sep2018.csv'
source('getData.R')

# look through each time slice and determine the distances from the data points
# to the mean in that time slice

dists <- NULL

for (tt in 1:n_time) {
    idx <- which(time[tt]-data_calib$age < 10 & time[tt]-data_calib$age >= 0)
    if (length(idx) > 1) {
        center <- mean(data_calib$co2[idx])
        dists <- c(dists, abs(data_calib$co2[idx]-center))
    }
}

print(mean(dists))
print(median(dists))
hist(log10(dists))


if(FALSE) {
# putting inverse gamma prior on variance, has units of (ppmv CO2)^2
# --> aim to center at median(dists^2) = 14055 ppmv^2

# some options:

# widest:
alph <- .1
beta <- 1/.000062
ftmp <- dinvgamma(x_co2^2, shape=alph, scale=1/beta); plot(x_co2, ftmp, col='red', type='l', xlim=c(0,3000))
x_co2[which.max(ftmp)]

# moderate width:
alph <- 1.1
beta <- 1/.000033
ftmp <- dinvgamma(x_co2^2, shape=alph, scale=1/beta); plot(x_co2, ftmp, col='red', type='l', xlim=c(0,3000))
x_co2[which.max(ftmp)]

# narrowest:
alph <- 2.1
beta <- 1/.000022
ftmp <- dinvgamma(x_co2^2, shape=alph, scale=1/beta); plot(x_co2, ftmp, col='red', type='l', xlim=c(0,3000))
x_co2[which.max(ftmp)]

# putting log-normal prior on stdev, has units of (ppmv CO2)
# --> aim to center at median(dists) = 119 ppmv and mean(dists) around 275 ppmv

logmean <- log(119) # know median = exp(logmean) for log-normal
logsd <- 1.3
ftmp <- dlnorm(x_co2, meanlog=logmean, sdlog=logsd); plot(x_co2, ftmp, col='red', type='l', xlim=c(0,3000))
print(exp(logmean+logsd*logsd/2))

}

##==============================================================================
## End
##==============================================================================
