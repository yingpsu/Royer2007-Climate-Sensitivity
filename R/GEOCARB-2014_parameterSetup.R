##==============================================================================
## GEOCARB-2014_parametersSetup.R
##
## questions? Tony Wong (twong@psu.edu)
##==============================================================================

l_Godderis <- TRUE  #set to "TRUE" to run time arrays of fA, fAw/fA, fD, and GEOG
                    #from Godderis et al, 2012; set to "FALSE" to run standard
                    #time arrays from GEOCARBSULF.
iteration_threshold <- 10 #maximum number of times the convergence equation for
                          #CO2 will iterate before signaling a failed run; in a
                          #test with reasonably-well-constrained input parameters
                          #(similar to simulations presented in Royer et al, 20XX),
                          #the number of iterations never exceeded 7; this variable
                          #is here mostly as a failsafe stop-gap.

# reading in the two input files

# input constant parameters and time-varying parameters
input <- read.csv("../input_data/GEOCARB_input_summaries_calib.csv")

# time-varying parameters
#y=young; a=old; p=pyrite; s=sulfate; c=carbonate; si=silicates;
# g=organic matter; b=burial; m=degassing; w=weathering
#masses are in units of 10^18 mol
#fluxes ("F" prefix) are in units of 10^18 mol Myrs-1
#rates ("k" prefix) are in units of Myrs-1
#stable isotopic compositions ("d" prefix) are in per mil units
time_arrays <- read.csv("../input_data/GEOCARB_input_arrays.csv")
ageN <- length(time_arrays$age) #number of time-steps
age <- matrix (time_arrays$age, nrow=ageN, ncol=1); colnames(age) <- "age (Myrs ago)"

# if using Godderis et al 2012 vlaues for fA, fAw_FA, fD, and GEOG...
if (l_Godderis) {
  arrays_Godderis <- c('fA','fAw_fA', 'fD', 'GEOG')
  for (i in 1:length(arrays_Godderis)) {
    # reset the array to match Godderis values
	time_arrays[,arrays_Godderis[i]] <- time_arrays[,paste(arrays_Godderis[i],'Godderis',sep='_')]
	# also reset the error associated with the array
	time_arrays[,paste('e',arrays_Godderis[i],sep='')] <- time_arrays[,paste(paste('e',arrays_Godderis[i],sep=''),'Godderis',sep='_')]
  }
}

ind_time <- which(input$type=='time array')
ind_const <- which(input$type=='constant')
parnames_time <- as.character(input$parameter[ind_time])
parnames_const <- as.character(input$parameter[ind_const])

ind_calib <- which(input$calib==1)
parnames_calib <- c(intersect(ind_calib, ind_time), intersect(ind_calib, ind_const))

# what are the columns within time_arrays of the calibration/fixed guys?
ind_time_calib <- intersect(ind_calib, ind_time)
ind_time_fixed <- setdiff(ind_time, ind_time_calib)

# what are the columns within input of the calibration/fixed guys?
ind_const_calib <- intersect(ind_calib, ind_const)
ind_const_fixed <- setdiff(ind_const, ind_const_calib)

# this needs to be a matrix initially, but for being passed around
# in the calibration, should be a vector (which matrix(..) can unwrap)
ind_time_mean <- NULL
ind_time_err <- NULL
for (i in 1:length(parnames_time)) {
  ind_time_mean <- c(ind_time_mean, match(parnames_time[i], colnames(time_arrays)))
  ind_time_err <- c(ind_time_err, match(paste('e',parnames_time[i],sep=''), colnames(time_arrays)))
}

par_time0 <- time_arrays[,ind_time_mean]
par_time_fixed0 <- data.frame(par_time0[,ind_time_fixed])
colnames(par_time_fixed0) <- parnames_time[ind_time_fixed]
par_time_calib0 <- data.frame(par_time0[,ind_time_calib])
colnames(par_time_calib0) <- parnames_time[ind_time_calib]

# unwrap par_time_fixed0 and calib0 into long vectors
par_time0_vec <- unlist(par_time0)
par_time_fixed0_vec <- unlist(par_time_fixed0)
par_time_calib0_vec <- unlist(par_time_calib0)

# unwrap the names corresponding to the values
parnames_time_fixed0_vec <- unlist(matrix(t(replicate(ageN, colnames(par_time_fixed0)))))
parnames_time_calib0_vec <- unlist(matrix(t(replicate(ageN, colnames(par_time_calib0)))))

# To wrap them back up as a matrix, just do (for example)
#par_time0 <- matrix(par_time0_vec, nrow=ageN)

par_const0 <- input$mean[ind_const]
par_const_fixed0 <- input$mean[ind_const_fixed]
par_const_calib0 <- input$mean[ind_const_calib]

par_calib0 <- c(par_const_calib0, as.numeric(par_time_calib0_vec))
par_fixed0 <- c(par_const_fixed0, as.numeric(par_time_fixed0_vec))

parnames_calib <- c( as.character(input$parameter[ind_const_calib]) ,
                     parnames_time_calib0_vec)
parnames_fixed <- c( as.character(input$parameter[ind_const_fixed]) ,
                     parnames_time_fixed0_vec )

# change the indices fed into the calibration to reflect which values within
# par_calib and par_fixed are constant parameters and which are time-varying
if (length(ind_const_calib)>0) {ind_const_calib <- 1:length(ind_const_calib)
} else                         {ind_const_calib <- NULL}
if (length(ind_const_fixed)>0) {ind_const_fixed <- 1:length(ind_const_fixed)
} else                         {ind_const_fixed <- NULL}

if (length(ind_time_calib)>0) {ind_time_calib <- (length(ind_const_calib)+1):length(par_calib0)
} else                        {ind_time_calib <- NULL}
if (length(ind_time_fixed)>0) {ind_time_fixed <- (length(ind_const_fixed)+1):length(par_fixed0)
} else                        {ind_time_fixed <- NULL}

##==============================================================================
## End
##==============================================================================
