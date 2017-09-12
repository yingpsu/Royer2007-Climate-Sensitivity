##==============================================================================
# This r code is written originally by Dana Royer for the publication Royer et
# al (20XX). The core GEOCARBSULF code (the second section of the code below) is
# translated mostly from the BASIC scripts written by Robert Berner in 2010 for
# the GEOCARBSULFvolc model (acquired by Dana Royer in October 2013), with updates
# from the FORTRAN scripts written by Jeffrey Park for Park & Royer (2011; see
# http://jparkcodes.blogspot.com/2011_06_01_archive.html); the most important
# component taken from the FORTRAN scripts is the convergence function for
# estimating CO2.
#
# For r novices, running the model is straightforward. First, download the r
# program if you don't already have it (http://www.r-project.org/). Next, set the
# working directory (e.g., desktop) (File->Change dir...); this is where all files
# will be read from and written to. Next, put the two input csv files in your working
# directory. Finally, copy the entire content of the r-code file (i.e., the file
# you are reading right now) and paste it into the r console. The summary csv file
# will be outputted to your working directory automatically. A full, standard run
# (10000 resamples) takes ~5 minutes to run on a PC with a 3.1 GHz processor and
# 8 GB of RAM.
#
# There are two input files: the default values reproduce the light gray uncertainty envelopes in Royer et al (20XX); to reproduce the dark gray envelopes, reduce all two sigma errors by 75%. "GEOCARB_input_summaries.csv" provides a summary of all input parameters (see also Appendices 1-3 in Royer et al, 20XX), including all necessary information for the time-invariant constants. "GEOCARB_input_arrays.csv" provides the time-step mean and two sigma values for the subset of parameters that are time arrays. Please note: if you open GEOCARB_input_summaries.csv in Excel, the "-inf" cells will probably be read in as an error; to correct, type <<'-inf>>.
# In the "input summaries" file:
#  parameters ending in "_Godderis" are only activated if "Godderis" is set to TRUE (third line in the code below); in this case, the corresponding non-Godderis parameter choices are ignored.
#  the "type" column differentiates time-invariant parameters (constants) from time-dependent parameters (arrays).
#  the "resample" column describes whether a particular parameter will undergo resampling for error analysis. If you wish to forgo all resampling, set resampleN to 1 (first line in the code below); this will override whatever is written in this particular column.
#  the "distribution" column describes whether the resampled distributions will follow a normal or lognormal distribution.
#  the "lower_limit" and "upper_limit" columns describe whether resamples should be clipped at particular threshold values; in these cases, a value just inside the limit (+/- 0.0001) is assigned; in the default model run, this rarely happens.
#
#Useful references:
#Berner RA. 2004. The Phanerozoic Carbon Cycle. Oxford University Press.
#Berner RA. 2006a. GEOCARBSULF: a combined model for Phanerozoic atmospheric O2 and CO2. Geochimica et Cosmochimica Acta 70: 5653-5664.
#Berner RA. 2006b. Inclusion of the weathering of volcanic rocks in the GEOCARBSULF model. American Journal of Science 306: 295-302.
#Berner RA. 2008. Addendum to "Inclusion of the weathering of volcanic rocks in the GEOCARBSULF model". American Journal of Science 308: 100-103.
#Godderis Y, Donnadieu Y, Lefebvre V, Le Hir G & Nardin E. 2012. Tectonic control of continental weathering, atmospheric CO2, and climate over Phanerozoic times. Comptes Rendus Geoscience 344: 652-662.
#Park J & Royer DL. 2011. Geologic constraints on the glacial amplification of Phanerozoic climate sensitivity. American Journal of Science 311: 1-26.
#Royer DL XXXXXXXXXXXXXXXXXXXXXXXXXX.
#
##==============================================================================
## 10 Mar 2017 -- Modified by Ying Cui and Tony Wong (twong@psu.edu)
##==============================================================================

  ######################################
  ### code for GEOCARBSULFvolc model ###
  ######################################
  GEOCARBSULFvolc_forMCMC <- function(Matrix_56, Matrix_12, age, ageN)
  {
	  #generate empty matrices for individual calculations of O2 and CO2
	  O2_resamples <- matrix(nrow=ageN, ncol=1)
	  CO2_resamples <- matrix(nrow=ageN, ncol=1)

      # match constant parameters to their Matrix_56 row using the row names
      # of Matrix_56
	  ACT <- Matrix_56['ACT']
	  ACTcarb <- Matrix_56['ACTcarb']
	  VNV <- Matrix_56['VNV']
	  NV <- Matrix_56['NV']
	  exp_NV <- Matrix_56['exp_NV']
	  LIFE <- Matrix_56['LIFE']
	  GYM <- Matrix_56['GYM']
	  FERT <- Matrix_56['FERT']
	  exp_fnBb <- Matrix_56['exp_fnBb']
	  deltaT2X <- Matrix_56['deltaT2X']
	  GLAC <- Matrix_56['GLAC']
	  J <- Matrix_56['J']
	  n <- Matrix_56['n']
	  Ws <- Matrix_56['Ws']
	  exp_fD <- Matrix_56['exp_fD']
	  Fwpa_0 <- Matrix_56['Fwpa_0']
	  Fwsa_0 <- Matrix_56['Fwsa_0']
	  Fwga_0 <- Matrix_56['Fwga_0']
	  Fwca_0 <- Matrix_56['Fwca_0']
	  Fmg_0 <- Matrix_56['Fmg_0']
	  Fmc_0 <- Matrix_56['Fmc_0']
	  Fmp_0 <- Matrix_56['Fmp_0']
	  Fms_0 <- Matrix_56['Fms_0']
	  Fwsi_0 <- Matrix_56['Fwsi_0']
	  Xvolc_0 <- Matrix_56['Xvolc_0']
	  CAPd13C_0 <- Matrix_56['CAPd13C_0']
	  CAPd34S_0 <- Matrix_56['CAPd34S_0']
	  oxygen_570 <- Matrix_56['oxygen_570']
	  Gy_570 <- Matrix_56['Gy_570']
	  Cy_570 <- Matrix_56['Cy_570']
	  Ca_570 <- Matrix_56['Ca_570']
	  Ssy_570 <- Matrix_56['Ssy_570']
	  Spy_570 <- Matrix_56['Spy_570']
	  dlsy_570 <- Matrix_56['dlsy_570']
	  dlcy_570 <- Matrix_56['dlcy_570']
	  dlpy_570 <- Matrix_56['dlpy_570']
	  dlpa_570 <- Matrix_56['dlpa_570']
	  dlgy_570 <- Matrix_56['dlgy_570']
	  dlga_570 <- Matrix_56['dlga_570']
	  Rcy_570 <- Matrix_56['Rcy_570']
	  Rca_570 <- Matrix_56['Rca_570']
	  Rv_570 <- Matrix_56['Rv_570']
	  Rg_570 <- Matrix_56['Rg_570']
	  Fob <- Matrix_56['Fob']
	  COC <- Matrix_56['COC']
	  Ga <- Matrix_56['Ga']
	  Ssa <- Matrix_56['Ssa']
	  Spa <- Matrix_56['Spa']
	  ST <- Matrix_56['ST']
	  dlst <- Matrix_56['dlst']
	  CT <- Matrix_56['CT']
	  dlct <- Matrix_56['dlct']
	  kwpy <- Matrix_56['kwpy']
	  kwsy <- Matrix_56['kwsy']
	  kwgy <- Matrix_56['kwgy']
	  kwcy <- Matrix_56['kwcy']

      # match time-varying paramteres to their Matrix_12 column using the column
      # names of Matrix 12
	  Sr <- Matrix_12[,'Sr']
	  d13C <- Matrix_12[,'d13C']
	  d34S <- Matrix_12[,'d34S']
	  fR <- Matrix_12[,'fR']
	  fL <- Matrix_12[,'fL']
	  fA <- Matrix_12[,'fA']
      fAw_fA <- Matrix_12[,'fAw_fA']
      fD <- Matrix_12[,'fD']
      GEOG <- Matrix_12[,'GEOG']
	  RT <- Matrix_12[,'RT']
	  fSR <- Matrix_12[,'fSR']
	  fC <- Matrix_12[,'fC']

	  Dt <- 10  #time-step (millions of years, Myrs)
	  oxygen_0 <- 38 #mass of atmospheric O2 at the present-day

      i <- 1    # no resampling - just a single iteration with a single parameter set
      RCO2 <- 10  #atmospheric CO2 (ratio between CO2 at time t to the Pleistocene mean [taken as 250 ppm]); this initial value is a place-holder (it is solved for explicitly)
	  #these variables are recalculated at each time-step
	  oxygen <- oxygen_570[i]
	  Gy <- Gy_570[i]
	  Cy <- Cy_570[i]
      Ca <- Ca_570[i]
	  Ssy <- Ssy_570[i]
	  Spy <- Spy_570[i]
	  dlsy <- dlsy_570[i]
	  dlcy <- dlcy_570[i]
	  dlpy <- dlpy_570[i]
	  dlpa <- dlpa_570[i]
	  dlgy <- dlgy_570[i]
	  dlga <- dlga_570[i]
	  dlsa <- (dlst[i]*ST[i]-(dlpy*Spy+dlsy*Ssy+dlpa*Spa[i]))/Ssa[i]  #d34S value of old CaSO4 sulfur
	  dlca <- (dlct[i]*CT[i]-(dlgy*Gy+dlcy*Cy+dlga*Ga[i]))/Ca  #d13C value of old crustal carbonate carbon
	  Rcy <- Rcy_570[i]
	  Rca <- Rca_570[i]

      # start the nested time loop (j)
	  for (j in 1:ageN) {
        failed_run <- FALSE  #flag for failed runs
		t <- age[j] #age of time-step, in Myrs ago

		# calculate factors that are influenced by glacial vs. non-glacial state ("GCM")
        # for glacial periods between 260 and 330 Myrs ago and between 35 and 0 Myrs ago;
        # BASIC code calls for a 270-340 Myrs ago interval, but see Fielding et al 2008
        # (GSA Special Paper 441: 343-354) for justification
		if ((t<=330 & t>=260) || t<=40) {
		  GCM <- GLAC[i]*deltaT2X[i]/log(2) #in GEOCARBSULF, deltaT2X = GCM*ln(2); called capital gamma in Berner (2004)
		} else {
		  GCM <- deltaT2X[i]/log(2)
        }

		# calculate factors related to vegetation type
		# vegetation = domination by non-vascular plants
        if (t<=570 & t>380) {
          # effect of plants on weathering rate at time (t) to the present-day
		  fE <- LIFE[i]
          # effect of CO2 on plant-assisted weathering for carbonates at time (t) to the present-day
		  fBB <- (1+ACTcarb[i]*GCM*log(RCO2)-ACTcarb[i]*Ws[i]*(t/570)+ACTcarb[i]*GEOG[j])*RCO2^exp_fnBb[i]
        }
		#vegetation = ramp-up to gymnosperm domination
        if (t<=380 & t>350) {
		  fE <- (GYM[i]-LIFE[i])*((380-t)/30)+LIFE[i]
		  fBB <- ((380-t)/30)*((1+ACTcarb[i]*GCM*log(RCO2)-ACTcarb[i]*Ws[i]*(t/570)+ACTcarb[i]*GEOG[j])*(2*RCO2/(1+RCO2))^FERT[i])+((t-350)/30)*(1+ACTcarb[i]*GCM*log(RCO2)-ACTcarb[i]*Ws[i]*(t/570)+ACTcarb[i]*GEOG[j])*RCO2^exp_fnBb[i]
        }
		if (t<=350)  {
		  fBB <- (1+ACTcarb[i]*GCM*log(RCO2)-ACTcarb[i]*Ws[i]*(t/570)+ACTcarb[i]*GEOG[j])*(2*RCO2/(1+RCO2))^FERT[i]
        }
		# vegetation = gymnosperm domination
        if (t<=350 & t>130)  {
		  fE <- GYM[i]
        }
		# vegetation = ramp-up to angiosperm domination
        if (t<=130 & t>80) {
		  fE <- (1-GYM[i])*((130-t)/50)+GYM[i]
        }
		#vegetation = angiosperm domination
        if (t<=80)  {
		  fE <- 1
        }

		# calculate source fluxes (weathering and degassing); see pp. 5660-5661
        # in Berner (2006a) for a discussion of the "Fwp", "Fws", and "Fwg" fluxes
		Fwpy <- fA[j]*fR[j]*kwpy[i]*Spy              #weathering flux of young pyrite
		Fwsy <- fA[j]*fD[j]*kwsy[i]*Ssy              #weathering flux of young CaSO4 sulfur
		Fwgy <- fA[j]*fR[j]*kwgy[i]*Gy               #weathering flux of young organic carbon
		Fwcy <- fA[j]*fD[j]*fL[j]*fE*fBB*kwcy[i]*Cy  #weathering flux of young carbonate carbon
		Fwpa <- fR[j]*Fwpa_0[i]                      #weathering flux of old pyrite
		Fwsa <- fA[j]*fD[j]*Fwsa_0[i]                #weathering flux of old CaSO4 sulfur
		Fwga <- fR[j]*Fwga_0[i]                      #weathering flux of old organic carbon
		Fwca <- fA[j]*fD[j]*fL[j]*fE*fBB*Fwca_0[i]   #weathering flux of old carbonate carbon
		Fmp <- fSR[j]*Fmp_0[i]                       #sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite
		Fms <- fSR[j]*Fms_0[i]                       #sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4
		Fmg <- fSR[j]*Fmg_0[i]                       #carbon degassing flux from volcanism, metamorphism, and diagenesis of organics
		Fmc <- fSR[j]*fC[j]*Fmc_0[i]                 #carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonate
		Fyop <- Fwpa+Fmp                             #degassing + weathering flux of pyrite
		Fyos <- Fwsa+Fms                             #degassing + weathering flux of CaSO4 sulfur
		Fyog <- Fwga+Fmg                             #degassing + weathering flux of organic carbon
		Fyoc <- Fwca+Fmc                             #degassing + weathering flux of carbonate carbon

		# calculate sink fluxes (burial)
		CAPd34S <- CAPd34S_0[i]*((oxygen/oxygen_0)^n[i]) #isotopic fractionation between sulfate sulfur and pyrite sulfur (see Berner 2006a)
		CAPd13C <- CAPd13C_0[i]+J[i]*(oxygen/oxygen_0-1) #isotopic fractionation between carbonate and organic matter (see Berner 2006a)
		Fbp <- (1/CAPd34S)*((d34S[j]-dlsy)*Fwsy+(d34S[j]-dlsa)*Fwsa+(d34S[j]-dlpy)*Fwpy+(d34S[j]-dlpa)*Fwpa+(d34S[j]-dlsa)*Fms+(d34S[j]-dlpa)*Fmp) #burial flux of pyrite
		Fbg <- (1/CAPd13C)*((d13C[j]-dlcy)*Fwcy+(d13C[j]-dlca)*Fwca+(d13C[j]-dlgy)*Fwgy+(d13C[j]-dlga)*Fwga+(d13C[j]-dlca)*Fmc+(d13C[j]-dlga)*Fmg) #burial flux of organic carbon
		Fbs <- Fwpy+Fwpa+Fwsy+Fwsa+Fms+Fmp-Fbp  #burial flux of CaSO4 sulfur
		Fbc <- Fwgy+Fwga+Fwcy+Fwca+Fmc+Fmg-Fbg  #burial flux of carbonate carbon

		#incorporate the volcanic/non-volcanic components (the "volc" in GEOCARBSULFvolc; see Berner 2006b & 2008)
		Roc <- (Sr[j]/10000)+0.7  #87Sr/86Sr of seawater as recorded in carbonates
		Rg <- Rg_570[i]-NV[i]*(1-fR[j]^(1/exp_NV[i]))  #87Sr/86Sr of non-volcanic silicates (same as "Rnv" in Berner 2006b)
		A <- ((Rv_570[i]-Roc)*fSR[j]*Fob[i])/(Fbc-Fwcy-Fwca) #note: this is always a negative number because Roc > Rv
		B <- Fwcy/(Fbc-Fwcy-Fwca)
		D <- Fwca/(Fbc-Fwcy-Fwca)
		E <- Fbc/(Fbc-Fwcy-Fwca)
		Xvolc <- (A+B*Rcy+D*Rca-E*Roc+Rg)/(Rg-Rv_570[i]) #fraction of total Ca and Mg silicate weathering drived from volcanic rocks
		fvolc <- (VNV[i]*Xvolc+1-Xvolc)/(VNV[i]*Xvolc_0[i]+1-Xvolc_0[i])  #volcanic weathering effect at time (t) relative to present-day

		# calculate oxygen mass for the time-step
		if (t<570)  {
		  oxygen <- oxygen+(Fbg+(15/8)*Fbp)*Dt-(Fwgy+Fwga+Fmg)*Dt-(15/8)*(Fwpy+Fwpa+Fmp)*Dt
        }

		# expression that summarizes the chemical weathering of silicates at time
        # (t) relative to the present-day in the absence of climatic effects (other
        # than the effect of CO2 on carbonate weathering); it is calculated by
        # dividing total silicate chemical weathering at time (t)--as determined
        # by mass balance between burial and degassing (numerator)--by the non-climatic
        # processes that affect silicate chemical weathering multiplied by the
        # present-day chemical weathering flux, Fwsi_0 (denominator) (see
        # equation 5.2 in Berner 2004) (called "fB" in BASIC scripts)
		fwsi_no_climate <- (Fbc-Fwcy-Fwca)/(((fAw_fA[j]*fA[j]*fD[j])^exp_fD[i])*fE*fR[j]*fvolc*Fwsi_0[i])

		# test for failed runs (negative numbers, NaN's, or oxygen <5% or >50%)
		if (min(B,D,E,fwsi_no_climate,Fwpy,Fwsy,Fwgy,Fwcy,Fwpa,Fwsa,Fwga,Fwca,
                Fmp,Fms,Fmg,Fmc,CAPd13C,CAPd34S,Fbp,Fbg,Fbs,Fbc,Spy,Ssy,Gy,Cy,
                Ca,Xvolc,fvolc)<=0 ||
            is.nan(fwsi_no_climate+A+B+D+E+Fwpy+Fwsy+Fwgy+Fwcy+Fwpa+Fwsa+Fwga+
                   Fwca+Fmp+Fms+Fmg+Fmc+CAPd13C+CAPd34S+Fbp+Fbg+Fbs+Fbc+Spy+Ssy+
                   Gy+Cy+Ca+dlpy+dlpa+dlsy+dlsa+dlgy+dlga+dlcy+dlca+Rcy+Rca+
                   Xvolc+fvolc+GCM+oxygen) || 100*(oxygen/(oxygen+143))<5 ||
            100*(oxygen/(oxygen+143))>50)  {failed_run <- TRUE}

		# calculate atmospheric CO2 (ppm) through iterative convergence (see
        # FORTRAN scripts of Park & Royer 2011 for details); this is done by
        # calculating the climatic factors that affect silicate weathering rates
        # normalized to the present-day (fwsi_climate, calculated below; called
        # "Fbbs" in BASIC scripts); fwsi_climate combines the expressions
        # "fBt(CO2)" and "fBb(CO2)" in Berner (2004) (see his equations 2.6-7,
        # 2.29, & 5.2); through inversion of fwsi_climate ("W", "V", & "X" below)
        # and comparison to fwsi_no_climate (calculated above), RCO2 can be determined;
        # see pp. 72-76 in Berner (2004) for details; iterative convergence is
        # necessary because RCO2--taken from the previous time-step--is needed
        # to calculate some of the dependent parameters; the initial calculation
        # of the new time-step for RCO2 is therefore not likely to be correct
		RCO2_old <- 2*RCO2  #this is here simply to ensure that the convergence isn't satisfied on the first step
		iteration_count <- 0
		if (t<=570 & t>380 & failed_run==FALSE)  {
		  while (abs(RCO2/RCO2_old-1) > 0.01) {
		    iteration_count <- iteration_count+1
			RCO2_old <- RCO2
			fwsi_climate <- (RCO2^(exp_fnBb[i]+ACT[i]*GCM))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			W <- ((exp_fnBb[i]+ACT[i]*GCM)*RCO2^(-exp_fnBb[i]+ACT[i]*GCM))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			V <- (RCO2^(exp_fnBb[i]+ACT[i]*GCM))*exp_fD[i]*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^-(1-exp_fD[i]))*(RT[j]*GCM/RCO2)*exp(-ACT[i]*Ws[i]*t/570)*exp(ACT[i]*GEOG[j])
			if (is.nan(fwsi_climate+W+V)==TRUE || iteration_count==iteration_threshold)  {
			  failed_run <- TRUE
			  break()
			}
			if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V))) {
			  RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V))  #damp the iteration to avoid overshoot
			} else {
			  RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
			}
		  } #end of while loop
		} #end of t>380 loop

        # the expressions for this time interval are more complex because the
        # effects of plants on weathering are linearly mixed in; this helps to
        # prevent model failure
		if (t<=380 & t>350 & failed_run==FALSE) {
		  while  (abs(RCO2/RCO2_old-1) > 0.01) {
		    iteration_count <- iteration_count+1
		    RCO2_old <- RCO2
			fwsi_climate_old <- (RCO2^(exp_fnBb[i]+ACT[i]*GCM))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			W_old <- (exp_fnBb[i]+ACT[i]*GCM)*RCO2^(-exp_fnBb[i]+ACT[i]*GCM)*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			V_old <- RCO2^(exp_fnBb[i]+ACT[i]*GCM)*exp_fD[i]*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^-(1-exp_fD[i]))*(RT[j]*GCM/RCO2)*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])

			fwsi_climate_new <- ((2^FERT[i])*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RCO2)^(-FERT[i]))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			W_new <- (2^FERT[i])*(FERT[i]+ACT[i]*GCM)*(RCO2^(FERT[i]+ACT[i]*GCM-1))*((1+RCO2)^-FERT[i])*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			V_new <- (-FERT[i]*(1+RCO2)^-(1+FERT[i]))*((2^FERT[i])*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[i])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			X_new <- exp_fD[i]*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^-(1-exp_fD[i]))*(RT[j]*GCM/RCO2)*(2^FERT[i]*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RCO2)^-FERT[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])

			fwsi_climate <- (t-350)/31*fwsi_climate_old+(381-t)/31*fwsi_climate_new
			Fw_v_x <- (t-350)/31*(W_old + V_old)+(381-t)/31*(W_new + V_new + X_new)

			if (is.nan(fwsi_climate_old + W_old + V_old + fwsi_climate_new + W_new + V_new + X_new + fwsi_climate + Fw_v_x)==TRUE  || iteration_count==iteration_threshold)  {
			  failed_run <- TRUE
			  break()
			}
			if (RCO2>((fwsi_climate - fwsi_no_climate)/Fw_v_x))  {
			  RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/Fw_v_x)  #damp the iteration to avoid overshoot
			} else {
			  RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
			}
		  } #end of while loop
		} #end of t<=380 & t>350 loop

		if (t<=350 & failed_run==FALSE) {
		  while  (abs(RCO2/RCO2_old-1) > 0.01) {
		    iteration_count <- iteration_count+1
			RCO2_old <- RCO2
			fwsi_climate <- ((2^FERT[i])*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RCO2)^(-FERT[i]))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			W <- (2^FERT[i])*(FERT[i]+ACT[i]*GCM)*(RCO2^(FERT[i]+ACT[i]*GCM-1))*((1+RCO2)^(-FERT[i]))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			V <- (-FERT[i]*(1+RCO2)^-(1+FERT[i]))*((2^FERT[i])*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^exp_fD[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			X <- exp_fD[i]*((1+RT[j]*GCM*log(RCO2)-RT[j]*Ws[i]*(t/570)+RT[j]*GEOG[j])^-(1-exp_fD[i]))*(RT[j]*GCM/RCO2)*(2^FERT[i]*RCO2^(FERT[i]+ACT[i]*GCM))*((1+RCO2)^-FERT[i])*exp(-ACT[i]*Ws[i]*(t/570))*exp(ACT[i]*GEOG[j])
			if (is.nan(fwsi_climate + W+V+X)==TRUE  || iteration_count==iteration_threshold)  {
			  failed_run <- TRUE
			  break()
			}
			if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V+X))) {
			  RCO2 <- RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V+X))  #damp the iteration to avoid overshoot
			} else {
			  RCO2 <- 0.2*RCO2  #damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
			}
		  } #end of while loop
		} #end of t<=350 loop

        # test for failed runs when CO2 < 150 ppm or > 50000 ppm
        if (RCO2<0.6 || RCO2>200) {
		  failed_run <- TRUE
        }

		# calculate new masses for the next time-step
		Spy <- Spy+(Fbp-Fwpy-Fyop)*Dt                    #mass of young pyrite sulfur
		Ssy <- Ssy+(Fbs-Fwsy-Fyos)*Dt                    #mass of young CaSO4 sulfur
		Gy <- Gy+(Fbg-Fwgy-Fyog)*Dt                      #mass of young crustal organic carbon
		Cy <- Cy+(Fbc-Fwcy-Fyoc)*Dt                      #mass of young crustal carbonate carbon
		Ca <- CT[i]-Gy-Ga[i]-Cy-COC[i]                   #mass of old crustal carbonate carbon

		# calculate new isotopic values for the next time-step
		dlpy <- dlpy+((d34S[j]-dlpy-CAPd34S)*Fbp/Spy)*Dt #d34S of young pyrite sulfur
		dlpa <- dlpa+(Fyop*(dlpy-dlpa)/Spa[i])*Dt        #d34S of old pyrite sulfur
		dlsy <- dlsy+((d34S[j]-dlsy)*Fbs/Ssy)*Dt         #d34S value of young CaSO4 sulfur
		dlsa <- dlsa+(Fyos*(dlsy-dlsa)/Ssa[i])*Dt        #d34S value of old CaSO4 sulfur
		dlgy <- dlgy+((d13C[j]-dlgy-CAPd13C)*Fbg/Gy)*Dt  #d13C of young organic matter
		dlga <- dlga+(Fyog*(dlgy-dlga)/Ga[i])*Dt         #d13C of old organic matter
		dlcy <- dlcy+((d13C[j]-dlcy)*Fbc/Cy)*Dt          #d13C value of young crustal carbonate carbon
		dlca <- dlca+(Fyoc*(dlcy-dlca)/Ca)*Dt            #d13C value of old crustal carbonate carbon
		Rcy <- Rcy+((Roc-Rcy)*Fbc/Cy)*Dt                 #87Sr/86Sr of young carbonates undergoing weathering over the time-step
		Rca <- Rca+((Rcy-Rca)*Fyoc/Ca)*Dt                #87Sr/86Sr of old carbonates undergoing weathering over the time-step

		# fill in the calculated values for O2 (%) and CO2 (ppm) (by default,
        # failed runs will be filled with "NAs")
		if (failed_run==FALSE) {
		  CO2_resamples[j] <- RCO2*250                   #RCO2 is multiplied by 250 (Pleistocene mean CO2) to convert to ppm units
		  O2_resamples[j] <- 100*(oxygen/(oxygen+143))   #converts O2 mass to O2 (%)
		}  else {
		  RCO2 <- 10                                     #reset RCO2 seed for next run
		}  #end of if..else loop

		# test for estimated oxygen at present-day to be between 19-23%, and
        # estimated CO2 at present-day to be between 200-300 ppm; if not, the
        # whole time-series is considered a failed run
		if (t==0 & (is.nan(oxygen+RCO2) || 100*(oxygen/(oxygen+143))<19 ||
           100*(oxygen/(oxygen+143))>23 || RCO2<0.8 || RCO2>1.2))  {
			CO2_resamples[,i] <- NA
			O2_resamples[,i] <- NA
		}   #end of if statement for t=0

      } #end of nested time loop (j)

      return (list(CO2=CO2_resamples, O2=O2_resamples))
}


##==============================================================================
## End
##==============================================================================




# below -- codes potentially useful for setting up the calibration (priors)



######################################################
### code for generating and filling input matrices ###
######################################################

#########
#function for generating resamples for time-dependent arrays, including the potential clipping of distributions
resamples_array <- function(x, name, row_position) {
  ###"Ying Cui" define ageN here because MCMC complains (two lines below)
  time_arrays <- read.csv("GEOCARB_input_arrays.csv")
  ageN <- length(time_arrays$age) #number of time-steps
  ###"Ying Cui" define resampleN <-1 here because MCMC complains (1 line below)
  resampleN <- 1
  ###"Ying Cui" define input here because MCMC complains
  input <- read.csv("GEOCARB_input_summaries_temp.csv")
  ###"Ying Cui" set Godderis TRUE
  l_Godderis <- TRUE  #set to "TRUE" to run time arrays of fA, fAw/fA, fD, and GEOG from Godderis et al, 2012; set to "FALSE" to run standard time arrays from GEOCARBSULF.

  x <- matrix(nrow=ageN, ncol=resampleN)  #create empty matrix
  col_num <- which(colnames(time_arrays)==name) #matches a column number in the input arrays file to the array in question

  #if resampling=1 or if resampling is turned off for the time-dependent array in question (a full matrix is generated but with resampleN mean values at each time-step)
  if (resampleN==1 || input[row_position, "resample"]==FALSE) {
    x <- matrix(time_arrays[, col_num], nrow=ageN, ncol=resampleN)  }

  #resampling following a normal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="gaussian" & resampleN>1) {
    for (i in 1:ageN) {
      temp_resample <- NULL
      temp_resample <- rnorm(n=resampleN, mean=time_arrays[i,col_num], sd=time_arrays[i,col_num+1]/2)
      temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
      temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
      x[i,] <- matrix(temp_resample, nrow=1, ncol=resampleN)
    } #end of time-sequence loop (i)
  }

  #resampling following a lognormal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="lognormal" & resampleN>1) {
    for (i in 1:ageN) {
      temp_resample <- NULL
      temp_resample <- rlnorm(n=resampleN, meanlog=log(time_arrays[i,col_num]), sdlog=log(time_arrays[i,col_num+1])/2)
      temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
      temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
      x[i,] <- matrix(temp_resample, nrow=1, ncol=resampleN)
    } #end of time-sequence loop (i)
  }
  x <- x
} #end of function "resamples_array"

#########
#function for generating resamples for constant parameters, including the potential clipping of distributions
resamples_constant <- function(x, row_position) {
  x <- matrix(nrow=1, ncol=resampleN)  #create empty matrix

  #if resampling=1 or if resampling is turned off for the constant parameter in question (a full matrix is filled with the mean value)
  if (resampleN==1 || input[row_position, "resample"]==FALSE)  {
    x <- matrix(input[row_position, "mean"], nrow=1, ncol=resampleN)  }

  #resampling following a normal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="gaussian" & resampleN>1) {
    temp_resample <- NULL
    temp_resample <- rnorm(n=resampleN, mean=input[row_position,"median"], sd=input[row_position,"two_sigma"]/2)
    temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
    temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
    x <- matrix(temp_resample, nrow=1, ncol=resampleN)
  }

  #resampling following a lognormal distribution
  if (input[row_position, "resample"]==TRUE & input[row_position, "distribution_type"]=="lognormal" & resampleN>1) {
    temp_resample <- NULL
    temp_resample <- rlnorm(n=resampleN, meanlog=log(input[row_position,"median"]), sdlog=log(input[row_position,"two_sigma"])/2)
    temp_resample <- replace(temp_resample,temp_resample<=input[row_position,"lower_limit"],input[row_position,"lower_limit"]+0.0001)
    temp_resample <- replace(temp_resample,temp_resample>=input[row_position,"upper_limit"],input[row_position,"upper_limit"]-0.0001)
    x <- matrix(temp_resample, nrow=1, ncol=resampleN)
  }
  x <- x
} #end of function "resamples_constant"
