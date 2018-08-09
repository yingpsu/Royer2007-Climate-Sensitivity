##==============================================================================
## GEOCARB-2014_processData_skewNormal_filter.R
##
## Read previous naive CO2 proxy data.
## Filter out anything that does not fit the Royer et al 2014 precalibration
## bounds of 150 <= CO2 <= 50000 ppmv
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

## Read the old data file with the skew-normal fits
filename.data <- '../input_data/CO2_Proxy_Foster2017_calib_SN-co2_06Jun2017.csv'
data_to_assim <- cbind( c("paleosols" , TRUE),
                        c("alkenones" , TRUE),
                        c("stomata"   , TRUE),
                        c("boron"     , TRUE),
                        c("liverworts", TRUE) )
source('GEOCARB-2014_getData.R')

## Find and filter those that do not pass the Royer et al 2014 laugh test
ind_filter <- which((data_calib$co2 < 150) | (data_calib$co2 > 50000))
data_calib <- data_calib[-ind_filter,]

## Write a filtered data file
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.out <- paste('../input_data/CO2_Proxy_Foster2017_calib_SN-co2-filtered_',today,'.csv', sep='')
write.csv(data_calib, file=filename.out)


##==============================================================================
## End
##==============================================================================
