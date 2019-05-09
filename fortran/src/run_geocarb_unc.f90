!================================================================================
!================================================================================
! run_geocarb: Fortran 90 file for GEOCARBSULFvolc model
! This file is written by Ying Cui (ying.cui@dartmouth.edu) and Tony Wong (twong@psu.edu)
!================================================================================
! This file is part of MCMC calibration for GEOCARBSULFvolc.
!================================================================================

!---------------------------------------------------------------------------------
subroutine run_geocarb_unc(Matrix_56, Matrix_12, age, ageN, iteration_threshold, CO2_out, O2_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Matrix_56  	56 constant parameters
! |     Matrix_12   12 time-varying parameters, each of length ageN
! |
! |    Parameters:
! |     age         age of time-step, in millions of years ago; 0-570 Ma
! |     ageN        number of time steps; each time step = 10; ageN=58
! |
! |
! | Outputs:
! |     CO2_out	    Modeled CO2 results [ppmv]
! |     O2_out	    Modeled O2 results [%]
!  =========================================================================

implicit none

integer, parameter :: i4b = Selected_Int_Kind(9)
integer, parameter :: i8b = Selected_Int_Kind(18)
integer, parameter :: SP = Selected_Real_Kind(6,37)
integer, parameter :: DP = Selected_Real_Kind(15,307)

! explicit input
integer, intent(IN) :: ageN                           ! time series length
integer, intent(IN) :: iteration_threshold            ! maximum number of iterations of convergence equation
real(DP), dimension(ageN), intent(IN) :: age          ! times of time series (myr ago)
real(DP), dimension(57), intent(IN) :: Matrix_56      ! 56 constant parameters to be evaluated
real(DP), dimension(58,ageN), intent(IN) :: Matrix_12 ! 12 time-series parameters to be evaluated, each of length ageN

! explicit output
real(DP), dimension(ageN), intent(OUT) :: CO2_out
real(DP), dimension(ageN), intent(OUT) :: O2_out

! local quantities
real(DP) :: var       ! variance statistical parameter (not used in here)
integer :: i          ! number of looping for time steps
real(DP) :: ACT       !Activation energy for dissolution of Mg- and Ca-silicates on land. 0.03 is equivalent to 20 kj/mol; 0.08 to 55 kj/mol; see Pg. 28 of Berner (2004)
real(DP) :: ACTcarb   !Activation energy for dissolution of carbonates on land
real(DP) :: VNV       !Rate of chemical weathering in volcanic silicate rocks relative to non-volcanic silicate rocks
real(DP) :: NV        !Coefficient relating physical erosion to the mean 87Sr/86Sr of non-volcanic silicate rocks
real(DP) :: exp_NV    !Exponent relating physical erosion to the mean 87Sr/86Sr of non-volcanic silicate rocks
real(DP) :: LIFE      !Rate of chemical weathering in a minimally-vegetated world relative to present-day (angiosperm-dominated world)
real(DP) :: GYM       !Rate of chemical weathering by gymnosperms relative to angiosperms
real(DP) :: FERT      !Exponent reflecting the fraction of vegetation whose growth is stimulated by elevated CO2; FERT is related to enhanced chemical weathering by the Michaelis-Menton expression [2RCO2/(1+RCO2)]^FERT
real(DP) :: exp_fnBb  !Exponent used to describe effect of climate on silicate or carbonate weathering in the absence of vascular plants at time (t) relative to the present-day
real(DP) :: deltaT2X  !Climate sensitivity; [compare this to G("GCM") factor, where deltaT2X=G*ln(2)]
real(DP) :: GLAC      !Factor by which deltaT2X changes during times with large continental icesheets
real(DP) :: J         !Coefficient used to calculate CAPd13C (called "alphac" in BASIC code), the stable carbon isotopic fractionation between shallow-marine carbonate and shallow marine organic matter; CAPd13C=27+J*(oxygen/38-1)
real(DP) :: n         !Exponent used to calculate CAPd34S (called "alphas" in BASIC code), the stable sulfur isotopic fractionation between shallow-marine CaSO4 sulfur and pyrite sulfur; CAPd34S = 35((oxygen/38)^n)
real(DP) :: Ws        !Effect on temperature from the linear increase in solar luminosity over time
real(DP) :: exp_fD    !Exponent that scales the dilution of dissolved HCO3- with runoff (fD)
real(DP) :: Fwpa_0    !Sulfate flux from oxidative weathering of old pyrite at present-day
real(DP) :: Fwsa_0    !Sulfate flux from weathering of CaSO4 sulfur at present-day
real(DP) :: Fwga_0    !Carbon flux from weathering of old sedimentary organic matter at present-day
real(DP) :: Fwca_0    !Carbon flux from weathering of old Ca and Mg carbonates at present-day
real(DP) :: Fmg_0     !Carbon degassing flux from volcanism, metamorphism, and diagenesis of organic matter at present-day
real(DP) :: Fmc_0     !Carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonates at present-day
real(DP) :: Fmp_0     !Sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite at present-day
real(DP) :: Fms_0     !Sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4 sulfur at present-day
real(DP) :: Fswi      !Weathering flux for all Ca and Mg silicates at present-day
real(DP) :: Xvolc_0   !Fraction of total Ca and Mg silicate weathering derived from volcanic rocks at present-day
real(DP) :: CAPd13C_0 !Stable carbon isotopic fractionation between shallow-marine carbonate and shallow-marine organic matter at present-day
real(DP) :: CAPd34S_0 !Stable sulfur isotopic fractionation between shallow-marine CaSO4 sulfur and pyrite sulfur at present-day
real(DP) :: oxygen_570 !Mass of atmospheric O2 (7.5 and 143 correspond to 5% and 50% O2) at 570 Myrs ago
real(DP) :: Gy_570    !Mass of young crustal organic carbon at 570 Myrs ago (value in Berner 2006a is a typo)
real(DP) :: Cy_570    !Mass of young crustal carbonate carbon at 570 Myrs ago (value in Berner 2006a is a typo)
real(DP) :: Ca_570    !Mass of old crustal carbonate carbon at 570 Myrs ago (value in Berner 2006a is a typo)
real(DP) :: Ssy_570   !Mass of young CaSO4 sulfur at 570 Myrs ago
real(DP) :: Spy_570   !Mass of young pyrite sulfur at 570 Myrs ago
real(DP) :: dlsy_570  !d34S of young CaSO4 sulfur at 570 Myrs ago
real(DP) :: dlcy_570  !d13C of young carbonate carbon at 570 Myrs ago
real(DP) :: dlpy_570  !d34S of young pyrite sulfur at 570 Myrs ago
real(DP) :: dlpa_570  !d34S of old pyrite sulfur at 570 Myrs ago
real(DP) :: dlgy_570  !d13C of young organic matter at 570 Myrs ago
real(DP) :: dlga_570  !d13C of old organic matter at 570 Myrs ago
real(DP) :: Rcy_570   !87Sr/86Sr of young carbonates undergoing weathering at 570 Myrs ago
real(DP) :: Rca_570   !87Sr/86Sr of old carbonates undergoing weathering at 570 Myrs ago
real(DP) :: Rv_570    !87Sr/86Sr of sub-aerial and submarine volcanic rocks at 570 Myrs ago
real(DP) :: Rg_570    !87Sr/86Sr of non-volcanic silicates at 570 Myrs ago
real(DP) :: Fob       !Ca and Mg flux between basalt and seawater
real(DP) :: COC       !Mass of carbon in ocean
real(DP) :: Ga        !Mass of old crustal organic carbon (value in Berner 2006a is a typo)
real(DP) :: Ssa       !Mass of old CaSO4 sulfur
real(DP) :: Spa       !Mass of old pyrite sulfur
real(DP) :: ST        !Mass of sulfur in oceans + "interacting rocks" (i.e., sulfur in rocks undergoing weathering, burial, etc.)
real(DP) :: dlst      !d34S of ST
real(DP) :: CT        !Mass of carbon in oceans + "interacting rocks" (i.e., carbon in rocks undergoing weathering, burial, etc.)
real(DP) :: dlct      !d13C of CT
real(DP) :: kwpy      !Rate constant expressing mass dependence for young pyrite sulfur
real(DP) :: kwsy      !Rate constant expressing mass dependence for young CaSO4 sulfur
real(DP) :: kwgy      !Rate constant expressing mass dependence for young organic matter weathering
real(DP) :: kwcy      !Rate constant expressing mass dependence for young carbonate weathering
real(DP) :: Roc

! declare local versions for the 12 time-varying parameters
real(DP), dimension(ageN) :: Sr        !Strontium isotopic composition of shallow-marine carbonate
real(DP), dimension(ageN) :: d13C      !Stable carbon isotopic composition of shallow-marine carbonate; called "DLCOC" in BASIC code
real(DP), dimension(ageN) :: d34S      !Stable sulfur isotopic composition of shallow-marine carbonate; called "DLSOC" in BASIC code
real(DP), dimension(ageN) :: fR        !Effect of relief on chemical weathering at time (t) relative to the present-day
real(DP), dimension(ageN) :: fL        !Land area covered by carbonates at time (t) relative to the present-day
real(DP), dimension(ageN) :: fA        !Land area at time (t) relative to the present-day
real(DP), dimension(ageN) :: fD        !Global river runoff at time (t) relative to the present-day in the absence of changes in solar luminosity and CO2
real(DP), dimension(ageN) :: fAw_fA    !Fraction of land area undergoing chemical weathering
real(DP), dimension(ageN) :: RT        !Coefficient of continental runoff versus temperature change, where fD=1+RT*(T-T0)
real(DP), dimension(ageN) :: GEOG      !Change in land mean surface temperature that is undergoing chemical weathering  at time (t) relative to the present-day in the absence of changes in solar luminosity and CO2 (see p. 35 in Berner 2004)
real(DP), dimension(ageN) :: fSR       !Seafloor creation rate at time (t) relative to the present-day
real(DP), dimension(ageN) :: fC        !Effect of carbonate content of subducting oceanic crust on CO2 degassing rate at time (t) relative to the present-day; 0.75 for >150 Ma, then ramps up linearly to 1 at present-day

! declare other variables in the model
real(DP) :: Dt        !time step (Myr)
real(DP) :: oxygen_0  !initial O2%
real(DP) :: Fbp
real(DP) :: Fbg
real(DP) :: Fwp
real(DP) :: Fwg
real(DP) :: oxygen
real(DP) :: Gy
real(DP) :: Cy
real(DP) :: Ca
real(DP) :: Ssy
real(DP) :: Spy
real(DP) :: dlsy
real(DP) :: dlcy
real(DP) :: dlpy
real(DP) :: dlpa
real(DP) :: dlgy
real(DP) :: lga
real(DP) :: dlsa
real(DP) :: dlca
real(DP) :: dlga
real(DP) :: Rcy
real(DP) :: Rca
real(DP) :: fE
real(DP) :: fBB
real(DP) :: A
real(DP) :: B
real(DP) :: D
real(DP) :: E
real(DP) :: W
real(DP) :: W_old
real(DP) :: W_new
real(DP) :: V
real(DP) :: V_old
real(DP) :: V_new
real(DP) :: X
real(DP) :: X_new
real(DP) :: S
real(DP) :: capd13C       !Stable carbon isotopic fractionation between shallow-marine carbonate and shallow-marine organic matter at any time t
real(DP) :: capd34S       !Stable sulfur isotopic fractionation between shallow-marine CaSO4 sulfur and pyrite sulfur at any time t
real(DP) :: t
real(DP) :: GCM
real(DP) :: Xvolc
real(DP) :: fbc
real(DP) :: fbs
real(DP) :: fmc
real(DP) :: fmg
real(DP) :: fmp
real(DP) :: fms
real(DP) :: fwsi_0
real(DP) :: fvolc
real(DP) :: Fw_v_x
real(DP) :: fwca
real(DP) :: fwcy
real(DP) :: fwga
real(DP) :: fwgy
real(DP) :: fwpa
real(DP) :: fwpy
real(DP) :: fwsa
real(DP) :: fwsi_climate
real(DP) :: fwsi_climate_new
real(DP) :: fwsi_climate_old
real(DP) :: fwsi_no_climate
real(DP) :: fwsy
real(DP) :: fyoc
real(DP) :: fyog
real(DP) :: fyop
real(DP) :: fyos
real(DP) :: RCO2
real(DP) :: RCO2_old
real(DP) :: Rg
integer :: iteration_count
LOGICAL :: failed_run

! match constant parameters to their Matrix_56 row using the row names of Matrix_56
! ... or *would* use the row names in the R version.
! this is tougher, hard-coded for now to the appropriate column, where
! the R wrapper must set up which is which appropriately.

ACT=Matrix_56(1)
ACTcarb=Matrix_56(2)
VNV=Matrix_56(3)
NV=Matrix_56(4)
exp_NV=Matrix_56(5)
LIFE=Matrix_56(6)
GYM=Matrix_56(7)
FERT=Matrix_56(8)
exp_fnBb=Matrix_56(9)
deltaT2X=Matrix_56(10)
GLAC=Matrix_56(11)
J=Matrix_56(12)
n=Matrix_56(13)
Ws=Matrix_56(14)
exp_fD=Matrix_56(15)
Fwpa_0=Matrix_56(16)
Fwsa_0=Matrix_56(17)
Fwga_0=Matrix_56(18)
Fwca_0=Matrix_56(19)
Fmg_0=Matrix_56(20)
Fmc_0=Matrix_56(21)
Fmp_0=Matrix_56(22)
Fms_0=Matrix_56(23)
Fwsi_0=Matrix_56(24)
Xvolc_0=Matrix_56(25)
CAPd13C_0=Matrix_56(26)
CAPd34S_0=Matrix_56(27)
oxygen_570=Matrix_56(28)
Gy_570=Matrix_56(29)
Cy_570=Matrix_56(30)
Ca_570=Matrix_56(31)
Ssy_570=Matrix_56(32)
Spy_570=Matrix_56(33)
dlsy_570=Matrix_56(34)
dlcy_570=Matrix_56(35)
dlpy_570=Matrix_56(36)
dlpa_570=Matrix_56(37)
dlgy_570=Matrix_56(38)
dlga_570=Matrix_56(39)
Rcy_570=Matrix_56(40)
Rca_570=Matrix_56(41)
Rv_570=Matrix_56(42)
Rg_570=Matrix_56(43)
Fob=Matrix_56(44)
COC=Matrix_56(45)
Ga=Matrix_56(46)
Ssa=Matrix_56(47)
Spa=Matrix_56(48)
ST=Matrix_56(49)
dlst=Matrix_56(50)
CT=Matrix_56(51)
dlct=Matrix_56(52)
kwpy=Matrix_56(53)
kwsy=Matrix_56(54)
kwgy=Matrix_56(55)
kwcy=Matrix_56(56)
var=Matrix_56(57)

! match time-varying paramteres to their Matrix_12 column using the column
! names of Matrix 12
Sr=Matrix_12(:,1)
d13C=Matrix_12(:,2)
d34S=Matrix_12(:,3)
fR=Matrix_12(:,4)
fL=Matrix_12(:,5)
fA=Matrix_12(:,6)
fAw_fA=Matrix_12(:,7)
fD=Matrix_12(:,8)
GEOG=Matrix_12(:,9)
RT=Matrix_12(:,10)
fSR=Matrix_12(:,11)
fC=Matrix_12(:,12)

Dt = 10.0       !time-step (millions of years, Myrs)
oxygen_0 = 38.0 !mass of atmospheric O2 at the present-day

RCO2 = 10.0     !atmospheric CO2 (ratio between CO2 at time t to the Pleistocene mean [taken as 250 ppm]); this initial value is a place-holder (it is solved for explicitly)

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

! ---------------------------------------------------------------------------------------------------------------------
! do loop starts
! start the nested time loop (i)

    do i=1,ageN

        failed_run = .FALSE.
        t = age(i) !age of time-step, in Myrs ago

        ! calculate factors that are influenced by glacial vs. non-glacial state ("GCM")
        ! for glacial periods between 260 and 330 Myrs ago and between 35 and 0 Myrs ago;
        ! BASIC code calls for a 270-340 Myrs ago interval, but see Fielding et al 2008
        ! (GSA Special Paper 441: 343-354) for justification
        if (((t .le. 330.0) .AND. (t .ge. 260.0)) .OR. (t .le. 40.0)) then
            GCM = GLAC*deltaT2X/log(2.0) !in GEOCARBSULF, deltaT2X = GCM*ln(2); called capital gamma in Berner (2004)
!!!            111 PRINT *, "GCM is: ", GCM
        else
            GCM = deltaT2X/log(2.0)
!!!            112 PRINT *, "GCM is: ", GCM
        end if

        ! calculate factors related to vegetation type
        ! vegetation = domination by non-vascular plants
        if (t .le. 570.0 .AND. t .gt. 380.0) then

          ! effect of plants on weathering rate at time (t) to the present-day
          fE = LIFE

          ! effect of CO2 on plant-assisted weathering for carbonates at time (t) to the present-day
          fBB = (1.0+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570.0)+ACTcarb*GEOG(i))*RCO2**exp_fnBb

        end if

        ! vegetation = ramp-up to gymnosperm domination
        if ((t .le. 380.0) .AND. (t .gt. 350.0)) then

          fE = (GYM-LIFE)*((380.0-t)/30.0)+LIFE

          ! GEOG is the change in avg land temp due to geography only, obtained via GCM runs
!!!          GEOG(i) = temp(i)-12.4
          fBB = ((380.0-t)/30.0)*((1.0+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570.0)+ACTcarb*GEOG(i))*(2.0*RCO2/(1.0+RCO2))**FERT) + &
                ((t-350.0)/30.0)*(1.0+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570.0)+ACTcarb*GEOG(i))*RCO2**exp_fnBb

        end if

        if (t .le. 350.0) then
          fBB = (1.0+ACTcarb*GCM*log(RCO2)-ACTcarb*Ws*(t/570.0)+ACTcarb*GEOG(i))*(2.0*RCO2/(1.0+RCO2))**FERT
        end if

        ! vegetation = gymnosperm domination
        if ((t .le. 350.0) .AND. (t .gt. 130.0))  then
          fE = GYM
!!!          113 PRINT *, "fE is: ", fE
        end if

        ! vegetation = ramp-up to angiosperm domination
        if ((t .le. 130.0) .AND. (t .gt. 80.0)) then
          fE = (1.0-GYM)*((130.0-t)/50.0)+GYM
        end if

        ! vegetation = angiosperm domination
        if (t .le. 80.0)  then
          fE = 1.0
        end if

        ! calculate source fluxes (weathering and degassing); see pp. 5660-5661
        ! in Berner (2006a) for a discussion of the "Fwp", "Fws", and "Fwg" fluxes
        Fwpy=fA(i)*fR(i)*kwpy*Spy              !weathering flux of young pyrite
        Fwsy=fA(i)*fD(i)*kwsy*Ssy              !weathering flux of young CaSO4 sulfur
        Fwgy=fA(i)*fR(i)*kwgy*Gy               !weathering flux of young organic carbon
        Fwcy=fA(i)*fD(i)*fL(i)*fE*fBB*kwcy*Cy  !weathering flux of young carbonate carbon
        Fwpa=fR(i)*Fwpa_0                      !weathering flux of old pyrite
        Fwsa=fA(i)*fD(i)*Fwsa_0                !weathering flux of old CaSO4 sulfur
        Fwga=fR(i)*Fwga_0                      !weathering flux of old organic carbon
        Fwca=fA(i)*fD(i)*fL(i)*fE*fBB*Fwca_0   !weathering flux of old carbonate carbon
        Fmp=fSR(i)*Fmp_0                       !sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite
        Fms=fSR(i)*Fms_0                       !sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4
        Fmg=fSR(i)*Fmg_0                       !carbon degassing flux from volcanism, metamorphism, and diagenesis of organics
        Fmc=fSR(i)*fC(i)*Fmc_0                 !carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonate
        Fyop=Fwpa+Fmp                          !degassing + weathering flux of pyrite
        Fyos=Fwsa+Fms                          !degassing + weathering flux of CaSO4 sulfur
        Fyog=Fwga+Fmg                          !degassing + weathering flux of organic carbon
        Fyoc=Fwca+Fmc                          !degassing + weathering flux of carbonate carbon

        ! calculate sink fluxes (burial)
        CAPd34S = CAPd34S_0*((oxygen/oxygen_0)**n)  !isotopic fractionation between sulfate sulfur and pyrite sulfur (see Berner 2006a)
        CAPd13C = CAPd13C_0+J*(oxygen/oxygen_0-1.0) !isotopic fractionation between carbonate and organic matter (see Berner 2006a)
        Fbp = (1.0/CAPd34S)*((d34S(i)-dlsy)*Fwsy+(d34S(i)-dlsa)*Fwsa + &
                (d34S(i)-dlpy)*Fwpy+(d34S(i)-dlpa)*Fwpa + &
                (d34S(i)-dlsa)*Fms+(d34S(i)-dlpa)*Fmp) !burial flux of pyrite
        Fbg = (1.0/CAPd13C)*((d13C(i)-dlcy)*Fwcy+(d13C(i)-dlca)*Fwca + &
                (d13C(i)-dlgy)*Fwgy+(d13C(i)-dlga)*Fwga + &
                (d13C(i)-dlca)*Fmc+(d13C(i)-dlga)*Fmg) !burial flux of organic carbon
        Fbs = Fwpy+Fwpa+Fwsy+Fwsa+Fms+Fmp-Fbp  !burial flux of CaSO4 sulfur
        Fbc = Fwgy+Fwga+Fwcy+Fwca+Fmc+Fmg-Fbg  !burial flux of carbonate carbon

        !incorporate the volcanic/non-volcanic components (the "volc" in GEOCARBSULFvolc; see Berner 2006b & 2008)
        Roc = (Sr(i)/10000.0)+0.7                     !87Sr/86Sr of seawater as recorded in carbonates
        Rg = Rg_570-NV*(1.0-fR(i)**(1.0/exp_NV))      !87Sr/86Sr of non-volcanic silicates (same as "Rnv" in Berner 2006b)
        A = ((Rv_570-Roc)*fSR(i)*Fob)/(Fbc-Fwcy-Fwca) !note: this is always a negative number because Roc > Rv
        B = Fwcy/(Fbc-Fwcy-Fwca)
        D = Fwca/(Fbc-Fwcy-Fwca)
        E = Fbc/(Fbc-Fwcy-Fwca)
        Xvolc = (A+B*Rcy+D*Rca-E*Roc+Rg)/(Rg-Rv_570)  !fraction of total Ca and Mg silicate weathering drived from volcanic rocks
        fvolc = (VNV*Xvolc+1.0-Xvolc)/(VNV*Xvolc_0+1.0-Xvolc_0)  !volcanic weathering effect at time (t) relative to present-day

        ! calculate oxygen mass for the time-step
        if (t .lt. 570.0)  then
          oxygen = oxygen+(Fbg+(15.0/8.0)*Fbp)*Dt-(Fwgy+Fwga+Fmg)*Dt-(15.0/8.0)*(Fwpy+Fwpa+Fmp)*Dt
        end if

        ! expression that summarizes the chemical weathering of silicates at time
        ! (t) relative to the present-day in the absence of climatic effects (other
        ! than the effect of CO2 on carbonate weathering); it is calculated by
        ! dividing total silicate chemical weathering at time (t)--as determined
        ! by mass balance between burial and degassing (numerator)--by the non-climatic
        ! processes that affect silicate chemical weathering multiplied by the
        ! present-day chemical weathering flux, Fwsi_0 (denominator) (see
        ! equation 5.2 in Berner 2004) (called "fB" in BASIC scripts)

        fwsi_no_climate=(Fbc-Fwcy-Fwca)/(((fAw_fA(i)*fA(i)*fD(i))**exp_fD)*fE*fR(i)*fvolc*Fwsi_0)

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

        RCO2_old = 2.0*RCO2  !this is here simply to ensure that the convergence isn't satisfied on the first step iteration_count=0
        iteration_count = 0

        if ( (t .le. 570.0) .AND. (t .gt. 380.0) .AND. (.not. failed_run))  then
          do while (abs((RCO2/RCO2_old)-1.0) .gt. 0.01)
            iteration_count = iteration_count+1
            RCO2_old = RCO2
            fwsi_climate = (RCO2**(exp_fnBb+ACT*GCM))*((1.0+RT(i)*GCM*log(RCO2) - &
                             RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            W = ((exp_fnBb+ACT*GCM)*RCO2**(-exp_fnBb+ACT*GCM))*((1.0+RT(i)*GCM*log(RCO2) - &
                  RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            V = (RCO2**(exp_fnBb+ACT*GCM))*exp_fD*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0) + &
                  RT(i)*GEOG(i))**(-(1.0-exp_fD)))*(RT(i)*GCM/RCO2)*exp(-ACT*Ws*t/570.0)*exp(ACT*GEOG(i))
            if (ISNAN(fwsi_climate+W+V) .OR. (iteration_count==iteration_threshold)) then
              failed_run=.TRUE.
              exit
            end if

            if (RCO2 .gt. ((fwsi_climate - fwsi_no_climate)/(W+V))) then
              RCO2 = RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V))  !damp the iteration to avoid overshoot
            else
              RCO2 = 0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
            end if
          end do ! end of while loop
        end if   ! end of t>380 loop

        ! the expressions for this time interval are more complex because the
        ! effects of plants on weathering are linearly mixed in; this helps to
        ! prevent model failure
        if ((t .le. 380.0) .AND. (t .gt. 350.0) .AND. (.not. failed_run)) then
          do while (abs((RCO2/RCO2_old)-1.0) .gt. 0.01)
            iteration_count = iteration_count+1
            RCO2_old = RCO2
            fwsi_climate_old = (RCO2**(exp_fnBb+ACT*GCM))*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0) + &
                                 RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            W_old = (exp_fnBb+ACT*GCM)*RCO2**(-exp_fnBb+ACT*GCM)*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0) + &
                      RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            V_old = RCO2**(exp_fnBb+ACT*GCM)*exp_fD*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0) + &
                      RT(i)*GEOG(i))**(-(1.0-exp_fD)))*(RT(i)*GCM/RCO2)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            fwsi_climate_new = ((2.0**FERT)*RCO2**(FERT+ACT*GCM))*((1.0+RCO2)**(-FERT))*((1.0+RT(i)*GCM*log(RCO2) - &
                                 RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            W_new = (2.0**FERT)*(FERT+ACT*GCM)*(RCO2**(FERT+ACT*GCM-1.0))*((1.0+RCO2)**(-FERT))*((1.0+RT(i)*GCM*log(RCO2) - &
                       RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            V_new = (-FERT*(1.0+RCO2)**(-(1.0+FERT)))*((2.0**FERT)*RCO2**(FERT+ACT*GCM))*((1.0+RT(i)*GCM*log(RCO2) - &
                       RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            X_new = exp_fD*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**(-(1.0-exp_fD))) * &
                       (RT(i)*GCM/RCO2)*(2.0**FERT*RCO2**(FERT+ACT*GCM))*((1.0+RCO2)**(-FERT)) * &
                       exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            fwsi_climate = (t-350.0)/31.0*fwsi_climate_old+(381.0-t)/31.0*fwsi_climate_new
            Fw_v_x = (t-350.0)/31.0*(W_old + V_old)+(381.0-t)/31.0*(W_new + V_new + X_new)

            if (ISNAN(fwsi_climate_old + W_old + V_old + fwsi_climate_new + W_new + V_new + X_new + fwsi_climate + Fw_v_x) .OR. &
                  iteration_count==iteration_threshold) then
              failed_run=.TRUE.
              exit
            end if

            if (RCO2 .gt. ((fwsi_climate - fwsi_no_climate)/Fw_v_x))  then
              RCO2=RCO2-0.9*((fwsi_climate - fwsi_no_climate)/Fw_v_x)  !damp the iteration to avoid overshoot
            else
              RCO2=0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
            end if
          end do  !end of while loop
        end if    !end of t<=380 & t>350 loop

        if ( (t .le. 350.0) .AND. (.not. failed_run) ) then
          do while( abs((RCO2/RCO2_old)-1.0) .gt. 0.01)
            iteration_count = iteration_count+1
            RCO2_old = RCO2
            fwsi_climate = ((2.0**FERT)*RCO2**(FERT+ACT*GCM))*((1.0+RCO2)**(-FERT))*((1.0+RT(i)*GCM*log(RCO2) - &
                             RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            W = (2.0**FERT)*(FERT+ACT*GCM)*(RCO2**(FERT+ACT*GCM-1.0))*((1.0+RCO2)**(-FERT)) * &
                   ((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            V = (-FERT*(1.0+RCO2)**(-(1.0+FERT)))*((2.0**FERT)*RCO2**(FERT+ACT*GCM)) * &
                   ((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**exp_fD)*exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            X = exp_fD*((1.0+RT(i)*GCM*log(RCO2)-RT(i)*Ws*(t/570.0)+RT(i)*GEOG(i))**(-(1.0-exp_fD))) * &
                   (RT(i)*GCM/RCO2)*(2.0**FERT*RCO2**(FERT+ACT*GCM))*((1.0+RCO2)**(-FERT)) * &
                   exp(-ACT*Ws*(t/570.0))*exp(ACT*GEOG(i))
            if (ISNAN(fwsi_climate + W+V+X) .OR. (iteration_count==iteration_threshold))  then
              failed_run=.TRUE.
              exit
            end if

            if (RCO2 .gt. ((fwsi_climate - fwsi_no_climate)/(W+V+X))) then
              RCO2 = RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V+X))  !damp the iteration to avoid overshoot
            else
              RCO2 = 0.2*RCO2  !damp the iteration (convert the iteration to geometric shrinkage to avoid nonpositive value in overshoot)
            end if

          end do !end of while loop
        end if   !end of t<=350 loop

        ! test for failed runs when CO2 < 150 ppm or > 50000 ppm
        if ( (RCO2 .lt. 0.6) .OR. (RCO2 .gt. 200.0) ) then
          failed_run = .TRUE.
        end if

        ! calculate new masses for the next time-step
        Spy=Spy+(Fbp-Fwpy-Fyop)*Dt                    !mass of young pyrite sulfur
        Ssy=Ssy+(Fbs-Fwsy-Fyos)*Dt                    !mass of young CaSO4 sulfur
        Gy=Gy+(Fbg-Fwgy-Fyog)*Dt                      !mass of young crustal organic carbon
        Cy=Cy+(Fbc-Fwcy-Fyoc)*Dt                      !mass of young crustal carbonate carbon
        Ca=CT-Gy-Ga-Cy-COC                            !mass of old crustal carbonate carbon

        ! calculate new isotopic values for the next time-step
        dlpy=dlpy+((d34S(i)-dlpy-CAPd34S)*Fbp/Spy)*Dt !d34S of young pyrite sulfur
        dlpa=dlpa+(Fyop*(dlpy-dlpa)/Spa)*Dt           !d34S of old pyrite sulfur
        dlsy=dlsy+((d34S(i)-dlsy)*Fbs/Ssy)*Dt         !d34S value of young CaSO4 sulfur
        dlsa=dlsa+(Fyos*(dlsy-dlsa)/Ssa)*Dt           !d34S value of old CaSO4 sulfur
        dlgy=dlgy+((d13C(i)-dlgy-CAPd13C)*Fbg/Gy)*Dt  !d13C of young organic matter
        dlga=dlga+(Fyog*(dlgy-dlga)/Ga)*Dt            !d13C of old organic matter
        dlcy=dlcy+((d13C(i)-dlcy)*Fbc/Cy)*Dt          !d13C value of young crustal carbonate carbon
        dlca=dlca+(Fyoc*(dlcy-dlca)/Ca)*Dt            !d13C value of old crustal carbonate carbon
        Rcy=Rcy+((Roc-Rcy)*Fbc/Cy)*Dt                 !87Sr/86Sr of young carbonates undergoing weathering over the time-step
        Rca=Rca+((Rcy-Rca)*Fyoc/Ca)*Dt                !87Sr/86Sr of old carbonates undergoing weathering over the time-step

        ! fill in the calculated values for O2 (%) and CO2 (ppm) (by default,
        ! failed runs will be filled with "NAs")
        CO2_out(i) = RCO2*250.0                       ! RCO2 is multiplied by 250 (Pleistocene mean CO2) to convert to ppm units
        O2_out(i) = 100.0*(oxygen/(oxygen+143.0))     ! converts O2 mass to O2 (%)
        if (failed_run) then
          !RCO2 = 10.0                                 ! reset RCO2 seed for next run
        end if

        ! test for estimated oxygen at present-day to be between 19-23%, and
        ! estimated CO2 at present-day to be between 200-300 ppm; if not, the
        ! whole time-series is considered a failed run
        if (((t==0.0) .AND. (ISNAN(oxygen+RCO2) .OR. (100.0*(oxygen/(oxygen+143.0)) .lt. 19.0) .OR. &
             (100.0*(oxygen/(oxygen+143.0)) .gt. 23.0) .OR. (RCO2 .lt. 0.8) .OR. (RCO2 .gt. 1.2))) .OR. &
              failed_run) then
          CO2_out(i) = 1.0/0.0
          O2_out(i) = 1.0/0.0
        end if

    end do  !end of nested time loop (i)

RETURN

END SUBROUTINE run_geocarb_unc
