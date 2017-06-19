!================================================================================
!================================================================================
! run_geocarb: Fortran 90 file for GEOCARBSULFvolc model
! This file is written by Ying Cui (ying.cui@dartmouth.edu) and Tony Wong (twong@psu.edu)
!================================================================================
! This file is part of MCMC calibration for GEOCARBSULFvolc.
!================================================================================

!---------------------------------------------------------------------------------
subroutine run_geocarb(Matrix_56, Matrix_12, age, ageN, CO2_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Matrix_56  	56 constant parameters
! |     Matrix_12   12 time-varying parameters
! |
! |    Parameters:
! |     age         age of time-step, in millions of years ago; 0-570 Ma
! |     ageN        number of time steps; each time step = 10; ageN=58
! |
! |
! | Outputs:
! |     CO2_output	Modeled CO2 results [ppmv]
!  =========================================================================


implicit none

integer(i4b), intent(IN) :: age, ageN ! time series length

! parameters
real(DP),     intent(IN) :: Matrix_56
real(DP),     intent(IN) :: Matrix_12


! generate empty matrices for individual calculations of O2 and CO2
! O2_resamples <- matrix(nrow=ageN, ncol=1)
! CO2_resamples <- matrix(nrow=ageN, ncol=1)

! match constant parameters to their Matrix_56 row using the row names
! of Matrix_56
Matrix_56="GEOCARB_input_summaries.txt"
open(8,Matrix_56,status="old") ! read in the constant 58 parameters
read(8,101) Matrix_56
Matrix_12="GEOCARB_input_arrays.txt"
open(8,"GEOCARB_input_arrays.txt",status="old") ! read in the constant 58 parameters
read(8,101) Matrix_12

ACT=Matrix_56('ACT')
ACTcarb=Matrix_56('ACTcarb')
VNV=Matrix_56('VNV')
NV=Matrix_56('NV')
exp_NV=Matrix_56('exp_NV')
LIFE=Matrix_56('LIFE')
GYM=Matrix_56('GYM')
FERT=Matrix_56('FERT')
exp_fnBb=Matrix_56('exp_fnBb')
deltaT2X=Matrix_56('deltaT2X')
GLAC=Matrix_56('GLAC')
J=Matrix_56('J')
n=Matrix_56('n')
Ws=Matrix_56('Ws')
exp_fD=Matrix_56('exp_fD')
Fwpa_0=Matrix_56('Fwpa_0')
Fwsa_0=Matrix_56('Fwsa_0')
Fwga_0=Matrix_56('Fwga_0')
Fwca_0=Matrix_56('Fwca_0')
Fmg_0=Matrix_56('Fmg_0')
Fmc_0=Matrix_56('Fmc_0')
Fmp_0=Matrix_56('Fmp_0')
Fms_0=Matrix_56('Fms_0')
Fwsi_0=Matrix_56('Fwsi_0')
Xvolc_0=Matrix_56('Xvolc_0')
CAPd13C_0=Matrix_56('CAPd13C_0')
CAPd34S_0=Matrix_56('CAPd34S_0')
oxygen_570=Matrix_56('oxygen_570')
Gy_570=Matrix_56('Gy_570')
Cy_570=Matrix_56('Cy_570')
Ca_570=Matrix_56('Ca_570')
Ssy_570=Matrix_56('Ssy_570')
Spy_570=Matrix_56('Spy_570')
dlsy_570=Matrix_56('dlsy_570')
dlcy_570=Matrix_56('dlcy_570')
dlpy_570=Matrix_56('dlpy_570')
dlpa_570=Matrix_56('dlpa_570')
dlgy_570=Matrix_56('dlgy_570')
dlga_570=Matrix_56('dlga_570')
Rcy_570=Matrix_56('Rcy_570')
Rca_570=Matrix_56('Rca_570')
Rv_570=Matrix_56('Rv_570')
Rg_570=Matrix_56('Rg_570')
Fob=Matrix_56('Fob')
COC=Matrix_56('COC')
Ga=Matrix_56('Ga')
Ssa=Matrix_56('Ssa')
Spa=Matrix_56('Spa')
ST=Matrix_56('ST')
dlst=Matrix_56('dlst')
CT=Matrix_56('CT')
dlct=Matrix_56('dlct')
kwpy=Matrix_56('kwpy')
kwsy=Matrix_56('kwsy')
kwgy=Matrix_56('kwgy')
kwcy=Matrix_56('kwcy')

! match time-varying paramteres to their Matrix_12 column using the column
! names of Matrix 12
Sr=Matrix_12('Sr')
d13C=Matrix_12('d13C')
d34S=Matrix_12('d34S')
fR=Matrix_12('fR')
fL=Matrix_12('fL')
fA=Matrix_12('fA')
fAw_fA=Matrix_12('fAw_fA')
fD=Matrix_12('fD')
GEOG=Matrix_12('GEOG')
RT=Matrix_12('RT')
fSR=Matrix_12('fSR')
fC=Matrix_12('fC')

Dt=10  !time-step (millions of years, Myrs)
oxygen_0=38 !mass of atmospheric O2 at the present-day

!i <- 1    # no resampling - just a single iteration with a single parameter set
RCO2=10  !atmospheric CO2 (ratio between CO2 at time t to the Pleistocene mean [taken as 250 ppm]); this initial value is a place-holder (it is solved for explicitly)
!these variables are recalculated at each time-step
oxygen = oxygen_570
Gy = Gy_570
Cy = Cy_570
Ca = Ca_570
Ssy = Ssy_570
Spy = Spy_570
dlsy = dlsy_570
dlcy = dlcy_570
dlpy = dlpy_570
dlpa = dlpa_570
dlgy = dlgy_570
dlga = dlga_570
dlsa = (dlst*ST-(dlpy*Spy+dlsy*Ssy+dlpa*Spa))/Ssa  !d34S value of old CaSO4 sulfur
dlca = (dlct*CT-(dlgy*Gy+dlcy*Cy+dlga*Ga))/Ca  !d13C value of old crustal carbonate carbon
Rcy = Rcy_570
Rca = Rca_570

! start the nested time loop (j)
    do j=1,ageN
!failed_run <- FALSE  #flag for failed runs
        t=age(j) !age of time-step, in Myrs ago

! calculate factors that are influenced by glacial vs. non-glacial state ("GCM")
! for glacial periods between 260 and 330 Myrs ago and between 35 and 0 Myrs ago;
! BASIC code calls for a 270-340 Myrs ago interval, but see Fielding et al 2008
! (GSA Special Paper 441: 343-354) for justification
    if ((t<=330 & t>=260) || t<=40) then
        GCM=GLAC*deltaT2X/log(2) !in GEOCARBSULF, deltaT2X = GCM*ln(2); called capital gamma in Berner (2004)
    else
        GCM=deltaT2X/log(2)
    end if

! calculate factors related to vegetation type
! vegetation = domination by non-vascular plants
    if (t<=570 & t>380) then
! effect of plants on weathering rate at time (t) to the present-day
    fE=LIFE
! effect of CO2 on plant-assisted weathering for carbonates at time (t) to the present-day
    fBB=(1+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570)+ACTcarb*GEOG)*RCO2^exp_fnBb
    end if

!vegetation = ramp-up to gymnosperm domination
    if (t<=380 & t>350) then
    fE=(GYM-LIFE)*((380-t)/30)+LIFE
    fBB=((380-t)/30)*((1+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570)+ACTcarb*GEOG(j))*(2*RCO2/(1+RCO2))^FERT)+((t-350)/30)*(1+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570)+ACTcarb*GEOG(j))*RCO2^exp_fnBb
    end if

    if (t<=350)  then
    fBB=(1+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570)+ACTcarb*GEOG(j))*(2*RCO2/(1+RCO2))^FERT
    end if

! vegetation = gymnosperm domination
    if (t<=350 & t>130)  then
    fE=GYM
    end if

! vegetation = ramp-up to angiosperm domination
    if (t<=130 & t>80) then
    fE=(1-GYM)*((130-t)/50)+GYM
    end if

! vegetation = angiosperm domination
    if (t<=80)  then
    fE=1
    end if

! calculate source fluxes (weathering and degassing); see pp. 5660-5661
! in Berner (2006a) for a discussion of the "Fwp", "Fws", and "Fwg" fluxes
Fwpy=fA(j)*fR(j)*kwpy*Spy              !weathering flux of young pyrite
Fwsy=fA(j)*fD(j)*kwsy*Ssy              !weathering flux of young CaSO4 sulfur
Fwgy=fA(j)*fR(j)*kwgy*Gy               !weathering flux of young organic carbon
Fwcy=fA(j)*fD(j)*fL(j)*fE*fBB*kwcy*Cy  !weathering flux of young carbonate carbon
Fwpa=fR(j)*Fwpa_0                      !weathering flux of old pyrite
Fwsa=fA(j)*fD(j)*Fwsa_0                !weathering flux of old CaSO4 sulfur
Fwga=fR(j)*Fwga_0                      !weathering flux of old organic carbon
Fwca=fA(j)*fD(j)*fL(j)*fE*fBB*Fwca_0   !weathering flux of old carbonate carbon
Fmp=fSR(j)*Fmp_0                       !sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite
Fms=fSR(j)*Fms_0                       !sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4
Fmg=fSR(j)*Fmg_0                       !carbon degassing flux from volcanism, metamorphism, and diagenesis of organics
Fmc=fSR(j)*fC(j)*Fmc_0                 !carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonate
Fyop=Fwpa+Fmp                             !degassing + weathering flux of pyrite
Fyos=Fwsa+Fms                             !degassing + weathering flux of CaSO4 sulfur
Fyog=Fwga+Fmg                             !degassing + weathering flux of organic carbon
Fyoc=Fwca+Fmc                             !degassing + weathering flux of carbonate carbon

! calculate sink fluxes (burial)
CAPd34S=CAPd34S_0*((oxygen/oxygen_0)^n) !isotopic fractionation between sulfate sulfur and pyrite sulfur (see Berner 2006a)
CAPd13C=CAPd13C_0+J*(oxygen/oxygen_0-1) !isotopic fractionation between carbonate and organic matter (see Berner 2006a)
Fbp=(1/CAPd34S)*((d34S(j)-dlsy)*Fwsy+(d34S(j)-dlsa)*Fwsa+(d34S(j)-dlpy)*Fwpy+(d34S(j)-dlpa)*Fwpa+(d34S(j)-dlsa)*Fms+(d34S(j)-dlpa)*Fmp) !burial flux of pyrite
Fbg=(1/CAPd13C)*((d13C(j)-dlcy)*Fwcy+(d13C(j)-dlca)*Fwca+(d13C(j)-dlgy)*Fwgy+(d13C(j)-dlga)*Fwga+(d13C(j)-dlca)*Fmc+(d13C(j)-dlga)*Fmg) !burial flux of organic carbon
Fbs=Fwpy+Fwpa+Fwsy+Fwsa+Fms+Fmp-Fbp  !burial flux of CaSO4 sulfur
Fbc=Fwgy+Fwga+Fwcy+Fwca+Fmc+Fmg-Fbg  !burial flux of carbonate carbon

!incorporate the volcanic/non-volcanic components (the "volc" in GEOCARBSULFvolc; see Berner 2006b & 2008)
Roc=(Sr(j)/10000)+0.7  !87Sr/86Sr of seawater as recorded in carbonates
Rg=Rg_570-NV*(1-fR(j)^(1/exp_NV))  !87Sr/86Sr of non-volcanic silicates (same as "Rnv" in Berner 2006b)
A=((Rv_570-Roc)*fSR(j)*Fob)/(Fbc-Fwcy-Fwca) !note: this is always a negative number because Roc > Rv
B=Fwcy/(Fbc-Fwcy-Fwca)
D=Fwca/(Fbc-Fwcy-Fwca)
E=Fbc/(Fbc-Fwcy-Fwca)
Xvolc=(A+B*Rcy+D*Rca-E*Roc+Rg)/(Rg-Rv_570) !fraction of total Ca and Mg silicate weathering drived from volcanic rocks
fvolc=(VNV*Xvolc+1-Xvolc)/(VNV*Xvolc_0+1-Xvolc_0)  !volcanic weathering effect at time (t) relative to present-day

! calculate oxygen mass for the time-step
    if (t<570)  then
    oxygen=oxygen+(Fbg+(15/8)*Fbp)*Dt-(Fwgy+Fwga+Fmg)*Dt-(15/8)*(Fwpy+Fwpa+Fmp)*Dt

! expression that summarizes the chemical weathering of silicates at time
! (t) relative to the present-day in the absence of climatic effects (other
! than the effect of CO2 on carbonate weathering); it is calculated by
! dividing total silicate chemical weathering at time (t)--as determined
! by mass balance between burial and degassing (numerator)--by the non-climatic
! processes that affect silicate chemical weathering multiplied by the
! present-day chemical weathering flux, Fwsi_0 (denominator) (see
! equation 5.2 in Berner 2004) (called "fB" in BASIC scripts)
    fwsi_no_climate=(Fbc-Fwcy-Fwca)/(((fAw_fA(j)*fA(j)*fD(j))^exp_fD)*fE*fR(j)*fvolc*Fwsi_0)
    end if

! test for failed runs (negative numbers, NaN's, or oxygen <5% or >50%)
    if (min(B,D,E,fwsi_no_climate,Fwpy,Fwsy,Fwgy,Fwcy,Fwpa,Fwsa,Fwga,Fwca, Fmp,Fms,Fmg,Fmc,CAPd13C,CAPd34S,Fbp,Fbg,Fbs,Fbc,Spy,Ssy,Gy,Cy,Ca,Xvolc,fvolc)<=0 || is.nan(fwsi_no_climate+A+B+D+E+Fwpy+Fwsy+Fwgy+Fwcy+Fwpa+Fwsa+Fwga+Fwca+Fmp+Fms+Fmg+Fmc+CAPd13C+CAPd34S+Fbp+Fbg+Fbs+Fbc+Spy+Ssy+Gy+Cy+Ca+dlpy+dlpa+dlsy+dlsa+dlgy+dlga+dlcy+dlca+Rcy+Rca+Xvolc+fvolc+GCM+oxygen) || 100*(oxygen/(oxygen+143))<5 ||100*(oxygen/(oxygen+143))>50) then
    failed_run=TRUE
    end if

! calculate atmospheric CO2 (ppm) through iterative convergence (see
! FORTRAN scripts of Park & Royer 2011 for details); this is done by
! calculating the climatic factors that affect silicate weathering rates
! normalized to the present-day (fwsi_climate, calculated below; called
! "Fbbs" in BASIC scripts); fwsi_climate combines the expressions
! "fBt(CO2)" and "fBb(CO2)" in Berner (2004) (see his equations 2.6-7,
! 2.29, & 5.2); through inversion of fwsi_climate ("W", "V", & "X" below)
! and comparison to fwsi_no_climate (calculated above), RCO2 can be determined;
! see pp. 72-76 in Berner (2004) for details; iterative convergence is
! necessary because RCO2--taken from the previous time-step--is needed
! to calculate some of the dependent parameters; the initial calculation
! of the new time-step for RCO2 is therefore not likely to be correct

RCO2_old=2*RCO2  !this is here simply to ensure that the convergence isn't satisfied on the first step
iteration_count=0
    end if !not sure if this end if statement is at the correct place

    if (t<=570 & t>380 & failed_run==FALSE)  then
        do while (abs(RCO2/RCO2_old-1) > 0.01)
        iteration_count=iteration_count+1
        RCO2_old=RCO2
        fwsi_climate=(RCO2^(exp_fnBb+ACT*GCM))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        W=((exp_fnBb+ACT*GCM)*RCO2^(-exp_fnBb+ACT*GCM))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        V=(RCO2^(exp_fnBb+ACT*GCM))*exp_fD*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^-(1-exp_fD))*(RT(j)*GCM/RCO2)*exp(-ACT*Ws*t/570)*exp(ACT*GEOG(j))

        if (is.nan(fwsi_climate+W+V)==TRUE || iteration_count==iteration_threshold)  then

        failed_run=TRUE
        break()
        end if

        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V))) then
        RCO2=RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V))  !damp the iteration to avoid overshoot
        else
        RCO2=0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        end if
        end do! end of while loop
    end if! end of t>380 loop

! the expressions for this time interval are more complex because the
! effects of plants on weathering are linearly mixed in; this helps to
! prevent model failure
    if (t<=380 & t>350 & failed_run==FALSE) then
        do while  (abs(RCO2/RCO2_old-1) > 0.01)
        iteration_count=iteration_count+1
        RCO2_old=RCO2
        fwsi_climate_old=(RCO2^(exp_fnBb+ACT*GCM))*((1+RT*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        W_old=(exp_fnBb+ACT*GCM)*RCO2^(-exp_fnBb+ACT*GCM)*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        V_old=RCO2^(exp_fnBb+ACT*GCM)*exp_fD*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^-(1-exp_fD))*(RT(j)*GCM/RCO2)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))

        fwsi_climate_new=((2^FERT)*RCO2^(FERT+ACT*GCM))*((1+RCO2)^(-FERT))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        W_new=(2^FERT)*(FERT+ACT*GCM)*(RCO2^(FERT+ACT*GCM-1))*((1+RCO2)^-FERT)*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        V_new=(-FERT*(1+RCO2)^-(1+FERT))*((2^FERT)*RCO2^(FERT+ACT*GCM))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG)^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        X_new=exp_fD*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^-(1-exp_fD))*(RT(j)*GCM/RCO2)*(2^FERT*RCO2^(FERT+ACT*GCM))*((1+RCO2)^-FERT)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))

        fwsi_climate=(t-350)/31*fwsi_climate_old+(381-t)/31*fwsi_climate_new
        Fw_v_x=(t-350)/31*(W_old + V_old)+(381-t)/31*(W_new + V_new + X_new)

        if (is.nan(fwsi_climate_old + W_old + V_old + fwsi_climate_new + W_new + V_new + X_new + fwsi_climate + Fw_v_x)==TRUE  || iteration_count==iteration_threshold) then
        failed_run=TRUE
        break()
        end if

        if (RCO2>((fwsi_climate - fwsi_no_climate)/Fw_v_x))  then
        RCO2=RCO2-0.9*((fwsi_climate - fwsi_no_climate)/Fw_v_x)  !damp the iteration to avoid overshoot
        else
        RCO2=0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        end if
        end do !end of while loop
    end if !end of t<=380 & t>350 loop

    if (t<=350 & failed_run==FALSE) then
        do while  (abs(RCO2/RCO2_old-1) > 0.01)
        iteration_count=iteration_count+1
        RCO2_old=RCO2
        fwsi_climate=((2^FERT)*RCO2^(FERT+ACT*GCM))*((1+RCO2)^(-FERT))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        W=(2^FERT)*(FERT+ACT*GCM)*(RCO2^(FERT+ACT*GCM-1))*((1+RCO2)^(-FERT))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        V=(-FERT*(1+RCO2)^-(1+FERT))*((2^FERT)*RCO2^(FERT+ACT*GCM))*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^exp_fD)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        X=exp_fD*((1+RT(j)*GCM*log(RCO2)-RT(j)*Ws*(t/570)+RT(j)*GEOG(j))^-(1-exp_fD))*(RT(j)*GCM/RCO2)*(2^FERT*RCO2^(FERT+ACT*GCM))*((1+RCO2)^-FERT)*exp(-ACT*Ws*(t/570))*exp(ACT*GEOG(j))
        if (is.nan(fwsi_climate + W+V+X)==TRUE  || iteration_count==iteration_threshold)  then
        failed_run=TRUE
        break()
        end if

        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V+X))) then
        RCO2=RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V+X))  !damp the iteration to avoid overshoot
        else
        RCO2=0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
        end if

        end do !end of while loop
    end if !end of t<=350 loop

! test for failed runs when CO2 < 150 ppm or > 50000 ppm
    if (RCO2<0.6 || RCO2>200) then
    failed_run=TRUE
    end if

! calculate new masses for the next time-step
Spy=Spy+(Fbp-Fwpy-Fyop)*Dt                    !mass of young pyrite sulfur
Ssy=Ssy+(Fbs-Fwsy-Fyos)*Dt                    !mass of young CaSO4 sulfur
Gy=Gy+(Fbg-Fwgy-Fyog)*Dt                      !mass of young crustal organic carbon
Cy=Cy+(Fbc-Fwcy-Fyoc)*Dt                      !mass of young crustal carbonate carbon
Ca=CT-Gy-Ga-Cy-COC                   !mass of old crustal carbonate carbon

! calculate new isotopic values for the next time-step
dlpy=dlpy+((d34S(j)-dlpy-CAPd34S)*Fbp/Spy)*Dt !d34S of young pyrite sulfur
dlpa=dlpa+(Fyop*(dlpy-dlpa)/Spa)*Dt        !d34S of old pyrite sulfur
dlsy=dlsy+((d34S(j)-dlsy)*Fbs/Ssy)*Dt         !d34S value of young CaSO4 sulfur
dlsa=dlsa+(Fyos*(dlsy-dlsa)/Ssa)*Dt        !d34S value of old CaSO4 sulfur
dlgy=dlgy+((d13C(j)-dlgy-CAPd13C)*Fbg/Gy)*Dt  !d13C of young organic matter
dlga=dlga+(Fyog*(dlgy-dlga)/Ga)*Dt         !d13C of old organic matter
dlcy=dlcy+((d13C(j)-dlcy)*Fbc/Cy)*Dt          !d13C value of young crustal carbonate carbon
dlca=dlca+(Fyoc*(dlcy-dlca)/Ca)*Dt            !d13C value of old crustal carbonate carbon
Rcy=Rcy+((Roc-Rcy)*Fbc/Cy)*Dt                 !87Sr/86Sr of young carbonates undergoing weathering over the time-step
Rca=Rca+((Rcy-Rca)*Fyoc/Ca)*Dt                !87Sr/86Sr of old carbonates undergoing weathering over the time-step

! fill in the calculated values for O2 (%) and CO2 (ppm) (by default,
! failed runs will be filled with "NAs")
!if (failed_run==FALSE) {
!CO2_resamples[j] <- RCO2*250                   #RCO2 is multiplied by 250 (Pleistocene mean CO2) to convert to ppm units
!O2_resamples[j] <- 100*(oxygen/(oxygen+143))   #converts O2 mass to O2 (%)
!}  else {
!RCO2 <- 10                                     #reset RCO2 seed for next run
!}  #end of if..else loop

! test for estimated oxygen at present-day to be between 19-23%, and
! estimated CO2 at present-day to be between 200-300 ppm; if not, the
! whole time-series is considered a failed run
if (t==0 & (is.nan(oxygen+RCO2) || 100*(oxygen/(oxygen+143))<19 ||
100*(oxygen/(oxygen+143))>23 || RCO2<0.8 || RCO2>1.2))  then
CO2 = NA
O2 = NA
end if

end !end of nested time loop (j)

return list(CO2=CO2_resamples, O2=O2_resamples)

end do
end if

RETURN

end subroutine run_geocarb
