c program gcsldata.f  
c fortran version of bob berner's GEOCARBSULF BASIC code 
c  this code developed in MacOS 10.4, using gcc compiler, 
c  takes 90minutes to run on a Macbook Pro with 2 GHz Intel Core Duo processor
c  -- looped over range of input parameters, set up to fit proxy data for CO2 ppm
c  11/22/06 JJP
c  proxy data has uncertainties that are sometimes large.  To make the data-fitting robust, 
c  compute chi-squared using an uncertainty in the log-domain.  
c  xf77 -o /Users/jjpark/bin/gcsldata_export gcsldata_export.f /Users/jjpark/Plotxy/plotlib.a
c
c   code uses plotit, a standalone fortran/C plotting subroutine in Xwindows
c  this code can be found on JPark's website http://earth.geology.yale.edu/~jjpark
c
c  code converted from BASIC to run grid search over parameters
c  Newton-Raphson damping applied here, to avoid NaN.
c  divergent transition over 380-350My is interpolated better, 
c  correcting earlier convergence error
c  the code reads a file proxydata.txt that contains proxy CO2 PPM data
c  the GEOCARBSULF carbon-cycle model is computed at time steps of 1-My from 570Ma to present
c  the proxy data ranges only from 420Ma to the present
c  
c  The heart of the code is five nested loops, from outermost to innermost
c  DeltaT = temp increase with CO2 doubling
c  ACT = dimensionless activation-energy for silicate weathering
c  FERT = plant CO2-fertilization coefficient
c  LIFE = liverwort-based weathering factor
c  GYM = gymnosperm-based weathering factor (1=angiosperm weathering parameter)
c  for each choice of parameters, we compute the carbon flux balance between terms 
c  that depend on silicate weathering (fBBS) and carbon burial/degassing (fB), 
c  each normalized with common factors chosen by Bob many years ago
c  The carbon-cycle balance depends on CO2 level, 
c  because the weathering of silicates accelerates with warmer temps induced by greenhouse gases
c
c  the code solves for the CO2 level that achieves the balance of carbon fluxes
c  The CO2-dependent parameter fB is calculated at each time step from carbon and
c  carbon isotope mass balance and values of all other parameters, based on
c  geological and biological data, that affect the carbon cycle. Then the value of
c  RCO2 (CO2 concentration normalized to the mean value for the past 1 million
c  years) is calculated by inversion from a complex expression for fB based on the
c  greenhouse effect, CO2 fertilization of plant weathering, solar evolution, and
c  changes in land temperature).
c  
c  c 
c  for detailed exposition of the theoretical GEOCARB models, see
c  Berner, R. A., The Phanerozoic Carbon Cycle: CO2 and O2, 150pp. Oxford University Press, 2004.
c  or a collection of GEOCARB and GEOCARBSULF papers published over the past decade by
c  Berner and Colleagues.
c
c  The code computes the aggregate data fit of GEOCARBSULF for each choice 
c  of the combined parameters, expressed as the fractional variance residual when expressed
c  in the log domain without uncertainty weighting.  
c  description of the proxy data can be found in 
c
c  Royer et al, CO2 as a primary driver of Phanerozoic climate, GSA Today, March 2004.
c
c  The data variance is expressed as the log-domain deviation from pre-industrial CO2-ppm 
c  the code stores all the variance residuals for all parameter combinations and computes the
c  percentiles of the model runs at each value of DeltaT.  The code also computes the fraction
c  of parameter space that achieves at most a chosen variance residual -- these fractions
c  are used to generate Bayesian PDFs.
c
c  percentile curves and empirical Bayesian PDFs are written in out1.dat and out2.dat
c
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      real*4 kwpy,kwgy,kwsy,kwcy,LIFE
      character*80 title,name   
      common/gcs1/DLSOC(600),DLCOC(600),alphac(600), GEOG(600)
      common/gcs2/fA(600),fSR(600),fD(600),fL(600),S(600)
      common/gcs3/temp(600),fAD(600),xx(600)
      common/gcs20/fA0(58),fSR0(72),fD0(58),fL0(58),S0(60) 
      common/gcs30/temp0(58),fAD0(60),DLS0(58),DLC0(58),al0(58)
      common/gcs4/Sr(600),R(600),DELTOT(600),DELRIV(600),Bas(600),
     x DELBAS(600),fR(600),SRBAS(600),sr0(68)
      common/loop/pp(10,10,10,10,60),dum(10000),weight(10000),wt(10)
      common/percento/p01(90),p10(90),p33(90),p50(90),
     x p66(90),p90(90),p99(90),ama(90)
      common/percentp/q01(90),q10(90),q33(90),q50(90),q66(90),q90(90),
     x q99(90),jkk(10,4,2)
      common/proxy/age(500),iage(500),pppm(500),dppm(500),res(500)
      common/dataf/gcsppm(600),aage(600)
      data fA0/1.0,1.0,1.0,1.0,0.91,0.92,0.94,0.87,0.79,0.81,
     x 0.82,0.85,0.88,0.86,0.84,0.83,0.85,0.87,0.89,0.91,
     x 0.91,0.90,0.90,0.88,0.89,0.86,0.84,0.81,0.81,0.80,
     x 0.80,0.80,0.79,0.77,0.76,0.73,0.70,0.69,0.69,0.72,
     x 0.74,0.74,0.74,0.62,0.63,0.67,0.62,0.62,0.62,0.65,
     x 0.65,0.66,0.64,0.66,0.85,0.90,0.95,1.00/
       data fD0/1.00,1.02,1.04,1.06,1.08,1.10,1.15,1.20,1.19,1.19,
     x 1.18,1.18,1.18,1.17,1.15,1.13,1.12,1.10,1.05,1.03,
     x 1.02,1.01,0.98,0.96,0.94,0.96,0.97,1.02,1.08,1.10,
     x 1.13,1.14,1.15,1.19,1.21,1.21,1.21,1.21,1.21,1.22,
     x 1.23,1.25,1.28,1.31,1.10,1.00,0.90,0.92,0.96,1.03,
     x 1.07,1.10,1.12,1.12,1.12,1.12,1.12,1.12/
       data temp0/12.4,12.4,12.4,12.5,12.5,12.5,11.6,11.0,11.2,11.3,
     x 11.5,11.7,11.9,12.0,12.4,12.8,13.1,13.4,13.7,14.0,
     x 14.2,14.2,14.1,14.1,14.1,14.3,14.2,13.9,13.7,13.5,
     x 13.3,13.1,13.0,13.0,12.9,13.0,13.1,13.3,13.4,13.5,
     x 12.5,11.5,10.5,12.5,15.0,17.0,18.2,17.0,16.3,15.5,
     x 14.5,13.5,13.0,13.0,13.0,12.0,12.0,12.0/
c  for the most recent 140My, the fSR parameter is given every 5My, not 10My
       data fSR0/1.00,0.98,0.98,0.98,0.98,1.02,1.04,1.11,1.22,1.41,
     x 1.46,1.41,1.35,1.30,1.35,1.37,1.33,1.30,1.37,1.52,
     x 1.65,1.76,1.74,1.67,1.63,1.67,1.52,1.26,1.27,1.27,
     x 1.27,1.21,1.14,1.09,1.10,1.14,1.12,1.10,1.09,1.11,
     x 1.14,1.10,1.07,1.21,1.21,1.14,1.07,1.07,1.28,1.43,
     x 1.43,1.40,1.38,1.36,1.38,1.48,1.53,1.55,1.55,1.52,
     x 1.52,1.50,1.49,1.50,1.52,1.62,1.70,1.59,1.35,1.20,1.04,1.00/
c FOR T= 0 TO 579
c IF T=<140 AND  INT(T/5)=T/5 THEN READ fSR(T)
c IF T>140 AND INT(T/10)=T/10 THEN READ fSR(T)
c IF T=<140 AND INT(T/5)<>T/5 THEN fSR(T)=fSR(T-1)
c IF T>140 AND INT(T/10)<>T/10 THEN fSR(T)=fSR(T-1)
c
c  data below are from Kampschulte and Strauss (2004)
       data DLS0/20.0,21.0,21.0,21.0,21.0,19.0,18.0,17.5,17.8,17.5,
     x 15.0,15.0,15.3,16.3,16.2,16.0,17.0,17.0,18.0,18.5,
     x 16.0,16.0,17.5,18.0,22.0,20.0,12.0,12.5,12.5,12.2,
     x 13.5,14.8,15.5,14.5,15.5,17.0,25.0,24.5,20.0,19.0,
     x 22.0,27.0,26.5,26.0,26.7,25.5,25.0,24.3,26.0,30.0,
     x 27.8,35.3,36.4,33.0,30.3,32.0,34.4,35.2/
c FOR T = 0 TO 579
c IF INT(T/10)=T/10 THEN READ DLSOC(T)
c IF INT (T/10)<>T/10 THEN DLSOC(T) =DLSOC(T-1)
c NEXT T      
c  following array IS smoothed fit TO maximum Veizer data with smoothed Silurian bump
c  also data of Korte (Tr),Popp, Mii (Perm+Carb),Krischvink(Camb), Saltzman (Ord, Sil), 
c  Lindh (Jurassic to 0)
      data DLC0/1.50,1.70,2.00,2.20,2.20,2.20,2.30,2.40,2.40,2.50,
     x 2.60,2.70,2.50,2.00,1.20,1.70,1.54,1.38,1.22,1.06,
     x 2.00,2.50,3.50,2.50,2.00,2.00,4.00,5.50,5.50,5.50,
     x 5.50,5.00,4.50,4.00,3.50,2.00,0.50,1.00,0.50,0.50,
     x 1.00,2.00,1.50,1.00,0.70,0.00,0.00,-1.0,-1.0,1.00,
     x -1.0,-0.5,0.00,0.00,-2.0,0.00,0.00,0.00/
c data are from Hayes ey al (1999)
      data al0/22.5,25.9,27.5,30.0,29.6,30.5,30.2,28.5,27.8,29.7,
     x 28.8,30.0,31.5,29.3,30.3,31.0,31.3,34.0,33.7,32.0,
     x 30.4,30.7,31.0,31.5,32.1,32.6,33.0,33.0,32.9,32.8,
     x 32.7,32.6,32.6,32.5,32.4,32.1,31.8,31.0,30.3,29.4,
     x 28.6,28.8,29.0,30.3,31.5,31.0,30.6,29.8,29.0,29.0,
     x 29.0,29.2,29.4,29.2,29.0,30.0,30.0,30.0/
c FOR T=  0 TO 579
c IF INT(T/10) =T/10 THEN READ alphac(T)
c IF INT(T/10)<>T/10 THEN  alphac(T)=alphac(T-1)
      data fL0/1.00,1.00,0.99,1.29,1.10,1.10,1.26,0.88,0.88,0.88,
     x 1.04,1.04,1.04,1.04,1.04,1.24,1.25,1.25,1.25,1.36,
     x 1.36,1.23,1.23,1.39,1.36,1.43,1.31,1.31,1.31,1.45,
     x 1.45,1.45,1.41,1.41,1.41,1.41,1.40,1.47,1.47,1.47,
     x 1.63,1.63,1.54,1.43,1.43,1.15,1.07,1.07,0.99,0.99,
     x 0.99,0.97,1.05,1.05,0.63,0.63,0.63,0.63/
c FOR T = 0 TO 579
c IF INT (T/10)=T/10 THEN READ fL(T)
c IF INT (T/10)<>T/10 THEN fL(T) = fL(T-1)
c sr0 is specified every 5my to 100Ma, prior times is 10My spacing
      data sr0/92.0,90.5,89.5,89.0,85.0,82.5,81.0,79.0,78.0,78.0,
     x 78.0,78.0,78.5,79.0,78.0,77.0,76.0,75.5,75.0,75.0,
     x 74.0,73.0,72.0,73.0,71.0,69.0,68.0,69.0,75.0,76.0,
     x 75.0,78.0,78.0,78.0,75.0,73.0,71.0,77.0,82.0,83.0,
     x 82.0,82.0,81.0,78.0,75.0,78.0,82.0,82.0,78.0,80.0,
     x 85.0,88.0,86.0,81.0,79.0,78.0,80.0,82.0,85.0,87.0,
     x 89.0,89.0,89.0,89.0,90.0,90.0,90.0,90.0/
c FOR T=0 TO 579
c IF T >100 AND INT(T/10)=T/10 THEN READ Sr(T)
c IF T>100 AND INT(T/10)<>T/10 THEN Sr(T)=Sr(T-1)
c IF T=<100 AND INT(T/5)=T/5 THEN READ Sr(T)
c IF T=<100 AND INT(T/5)<>T/5 THEN Sr(T)=Sr(T-1)


c  read proxydata file
c  uncertainties provided by Dana Royer are in a range, a max and min.  
c  I compute a factor-of-X uncertainty by dividing the upper bound by the value
c  and convert the 1-sigma into logarithmic misfit
      name='proxydata.txt'
      ndat=488
      open(8,file=name,form='formatted')
      read(8,101) name
  101 format(a)
      do i=1,ndat
        read(8,*) age(i),pppm(i),plo,phi
	dppm(i)=alog(phi/pppm(i))
	dum(i)=alog(pppm(i))
	iage(i)=1+(age(i)+0.5001)
      end do
      close(8)
      blob=0.
c      call plotit(age,dum,dppm,ndat,'proxydata','age(Ma)',
c     x 'log CO\\sub{2} ppm',
c     x  3,0,0.1,1,21)
c      call plotit_axes(0.,75.,0.,0.)
c      call plotit(age,dum,dppm,ndat,'proxydata','log age(Ma)',
c     x 'log CO\\sub{2} ppm',
c     x  3,0,0.1,1,22)
c      pause
c      call plotit_axes(0.,0.,0.,0.)
c      if(blob.eq.1.0) then
c  compute 10Ma moving average of CO2 PPM estimates
c  the average is weighted by inverse-variance in the log domain.
c It is interesting that the PPM estimates since the Eocene are so much
c lower than GEOCARB.  If the proxy data are unbiased, this implies that
c  there are more factors that control CO2 than GEOCARB includes -- an
c inference that we probably would not dismiss.  Short-term disagreement
c doesnt invalidate the argument about the longer-term changes.
      jdat=0
      do i=1,42
        sum=0.
	nsum=0
	dsum=0
	do j=1,ndat
	  if((iage(j)-1)/10+1.eq.i) then
	    nsum=nsum+1
	    sum=sum+alog(pppm(j))/(dppm(j)**2)
	    dsum=dsum+1.0/(dppm(j)**2)
	  endif
	end do
	if(nsum.gt.0) then
	  jdat=jdat+1
  	  dum(jdat)=exp(sum/dsum)
	  dum(42+jdat)=sqrt(1.0/dsum)
	  dum(jdat+126)=10.*i-4.
	endif
c  compute the chi-square of the moving average in log domain
c  can scale up the variance of the moving-average for large chi-square
        sum=0.
	do j=1,ndat
	  if((iage(j)-1)/10+1.eq.i) then
            sum=sum+(alog(dum(jdat)/pppm(j))/dppm(j))**2
	  endif
	end do
	if(nsum.gt.0) then
          sum=sum/nsum
          dum(84+jdat)=sum
	  if(nsum.eq.1) dum(84+jdat)=sum+1.
c	  print *,dum(jdat),dum(42+jdat),sum
	  if(sum.gt.1.0) dum(42+jdat)=dum(42+jdat)*sqrt(sum)
	endif
c  INSTEAD: compute the sample variance of the PPM estimates within the moving average 
c  in log domain -- dont use uncertainty weighting, so we are using log-ppm units
        sum=0.
	do j=1,ndat
	  if((iage(j)-1)/10+1.eq.i) then
            sum=sum+(alog(dum(jdat)/pppm(j)))**2
	  endif
	end do
	if(nsum.gt.1) then
          sum=sum/(nsum-1)
          dum(42+jdat)=sum
	  print *,dum(jdat),dum(42+jdat),sum
	endif
      end do
c      call plotit(age,pppm,dm,ndat,'proxydata','age(Ma)',
c     x 'CO\\sub{2} ppm',
c     x  2,2,0.1,1,0)
      ndat=jdat
      do i=1,ndat
        iage(i)=i*10-4
	pppm(i)=dum(i)
	dppm(i)=dum(42+i)
	age(i)=dum(126+i)
	dum(i)=pppm(i)*exp(dppm(i))
	dum(i+42)=pppm(i)*exp(-dppm(i))
c	print *,iage(i),pppm(i),dppm(i),age(i)
      end do
c      call plotit(age,dum,dm,ndat,'proxydata','age(Ma)',
c     x 'CO\\sub{2} ppm',
c     x  2,2,0.0,0,0)
c      call plotit(age,dum(43),dm,ndat,'proxydata','age(Ma)',
c     x 'CO\\sub{2} ppm',
c     x  2,2,0.0,0,0)  
c      call plotit(age,pppm,dm,ndat,'proxydata','age(Ma)',
c     x 'CO\\sub{2} ppm',
c     x  2,2,0.0,0,21)
c      call plotit_axes(0.,0.,0.,0.)
c      call plotit(age,dum(85),dm,ndat,
c     x 'average variance/dof in moving avg','age(Ma)',
c     x 'CO\\sub{2} ppm',2,2,0.1,1,22)
c      endif
c GEOCARB III + Rapid  Recyc, Hayes alphac, n for alphas, 
c NO MIMI. Korte Triassic
c  McDougall 87Sr at P/YTr boundary
c vary delta T from 0.4 to 10.1 for each value of ACT, FERT,LIFE and GYM
      ndel=82
      do iiit=1,ndel
      deltaT = 0.28+iiit*0.12
      GCM = deltaT/0.69
      varredmin=1000.
      do k1=10,1,-1
        ACT=0.03+(k1-1)*0.10/9.
      do k2=10,1,-1
c  limit plant response to CO2 enrichment to range of 20% to 80% efficiency
        FERT=0.2+(k2-1)*0.2/3.
      do k3=10,1,-1
        LIFE=0.125+(k3-1)*0.375/9.
      do k4=10,1,-1
        GYM=0.5+(k4-1)*0.7/9.
      Dt=1
      Ws=7.4
      SOC=38
      COC=2
      O=25.
      RCO2 =10.
      kwpy=0.01							
      kwsy=0.01							
      Spy=20.
      Spa=280
      Ssy=150
      Ssa=150
      kwgy=0.018							
      kwcy=0.018							
      Gy=250
      Ga=1000
      Cy=1000
      Ca=4000
      Fwpa1=0.25
      Fwsa1=0.5
      Fwga1=0.5
      Fwca1=2
      Fmg1=1.25
      Fmc1=6.67
      Fmp1=0.25
      Fms1=0.5
      St=600
      dlst=4
      CT=6252
      dlct=-3.5
      dlsy=35
      dlcy=3
      dlpy=-10
      dlpa=-10
      dlsa= (dlst*St-(dlpy*Spy+dlsy*Ssy+dlpa*Spa))/Ssa
      dlgy= -23.5
      dlga=-23.5
      dlca= (dlct*CT-(dlgy*Gy  +dlcy*Cy  +dlga*Ga))/Ca
      ZM=12.5							
      RBAS=.703
      FBASO=.92
      RRIV=.711
      Dt=1
      FRIV=3.37
c  fill in the Sr array
      do i=1,20
        jj=1+(i-1)*5
	do j=0,4
	  Sr(jj+j)=(j*sr0(i+1)+(5-j)*sr0(i))/5.
	end do
      end do
      do i=21,67
        jj=101+(i-21)*10
	do j=0,9
	  Sr(jj+j)=(j*sr0(i+1)+(10-j)*sr0(i))/10.
	end do
      end do
      Sr(571)=sr0(68)
c  fill in the fSR array
      do i=1,28
        jj=1+(i-1)*5
	do j=0,4
	  fSR(jj+j)=(j*fSR0(i+1)+(5-j)*fSR0(i))/5.
	end do
      end do
      do i=29,71
        jj=141+(i-29)*10
	do j=0,9
	  fSR(jj+j)=(j*fSR0(i+1)+(10-j)*fSR0(i))/10.
	end do
      end do
      fSR(571)=fSR0(72)
      do it=571,1,-1
c test case would be FSR(T)=1
        FBAS=FBASO*fSR(it)
        IF(it.lt.571) Bas(it)=Bas(it+1)+DELTOT(it+1)
        IF(it.eq.571) Bas(it)=0.709
        IF(it.eq.571) fR(it)=1
        DELRIV(it)=((RRIV-Bas(it))/ZM)*FRIV
        DELBAS(it)=((RBAS-Bas(it))/ZM)*FBAS
        DELTOT(it)=DELBAS(it)+DELRIV(it)
        R(it)=0.7+Sr(it)/10000.
        SRBAS(it)=(Bas(it)-0.7)*10000
        L=2
        fR(it)=1-L*(1-(Sr(it)/SRBAS(it)))
      end do
c  fill in the other arrays
      do i=1,57
        jj=1+(i-1)*10
	do j=0,9
	  fA(jj+j)=(j*fA0(i+1)+(10-j)*fA0(i))/10.
	  fD(jj+j)=(j*fD0(i+1)+(10-j)*fD0(i))/10.
	  temp(jj+j)=(j*temp0(i+1)+(10-j)*temp0(i))/10.
	  DLSOC(jj+j)=(j*DLS0(i+1)+(10-j)*DLS0(i))/10.
	  DLCOC(jj+j)=(j*DLC0(i+1)+(10-j)*DLC0(i))/10.
	  alphac(jj+j)=(j*al0(i+1)+(10-j)*al0(i))/10.
	  fL(jj+j)=(j*fL0(i+1)+(10-j)*fL0(i))/10.
	end do
      end do
      fA(571)=fA0(58)
      fD(571)=fD0(58)
      temp(571)=temp0(58)
      DLSOC(571)=DLS0(58)
      DLCOC(571)=DLC0(58)
      alphac(571)=al0(58)
      fL(571)=fL0(58)
c  OK, we finally have all initialized
      do i=1,579
        xx(i)=1-i
      end do
      blob=0.
c      open(8,file='out.out',form='formatted')
c        write(8,*) GYM,LIFE,ACT,FERT
      do i=571,1,-1
	fac=(i-1)/570.
        if(i.gt.341) then
	  RT=0.025
        elseif(i.le.341.and.i.gt.261) then 
	  RT=0.045
        elseif(i.le.261.and.i.gt.41) then 
          RT=0.025
        else
	  RT=0.045
	endif
        GAS =0.75
        if(i.gt.151) then 
          Fc=GAS
        else
	  Fc=(GAS)+((1.0-GAS)/150.)*(151-i)
	endif  
        if(i.gt.381) then  
	  fE=LIFE
        elseif(i.le.381.and.i.gt.351) then  
	  fE=(GYM-LIFE)*((381.-i)/30.)+LIFE
        elseif(i.le.351.and.i.gt.131) then   
	  fE=GYM
        elseif(i.le.131.and.i.gt.81) then  
	  fE=(1.-GYM)*((131-i)/50.)+GYM
        else
	  fE=1
	endif
        GEOG(i) = temp(i)-12.4
        if(i.eq.571) then
	  fBB=1
        elseif(i.le.570.and.i.gt.381) then  
          fBB=(1+0.087*GCM*alog(RCO2)-0.087*Ws*fac
     x                        +0.087*GEOG(i))*sqrt(RCO2)
        
        elseif(i.le.381.and.i.gt.351) then  
           fBB=((381.-i)/30.)*((1.0+0.087*GCM*alog(RCO2)-0.087*Ws*fac
     x  +0.087*GEOG(i))*(2.0*RCO2/(1.0+RCO2))**FERT)
     x  +((i-351)/30.)*(1.0+0.087*GCM*alog(RCO2)
     x  -0.087*Ws*fac+0.087*GEOG(i))*sqrt(RCO2)
        else 
c	  print *,fBB,GCM,RCO2,Ws,fac,FERT,GEOG(i)
          fBB=(1.0+0.087*GCM*alog(RCO2)-0.087*Ws*fac
     x 	    +0.087*GEOG(i))*(2.0*RCO2/(1.0+RCO2))**FERT
	endif
	Fwpy=fA(i)*fR(i)*kwpy*Spy
	Fwsy=fA(i)*fD(i)*kwsy*Ssy
	Fwgy=fA(i)*fR(i)*kwgy*Gy
	Fwcy=fA(i)*fD(i)*fL(i)*fE*fBB*kwcy*Cy
	Fwpa=fR(i)*Fwpa1
	Fwsa=fA(i)*fD(i)*Fwsa1
	Fwga=fR(i)*Fwga1
	Fwca=fA(i)*fD(i)*fL(i)*fE*fBB*Fwca1
	Fmp=fSR(i)*Fmp1
	Fms=fSR(i)*Fms1
	Fmg=fSR(i)*Fmg1
	Fmc=fSR(i)*Fc*Fmc1
	Fyop=Fwpa+Fmp
	Fyos=Fwsa+Fms
	Fyog=Fwga+Fmg
	Fyoc=Fwca+Fmc
	zzn=1.5
	J=3
	alphas =35.0*((O/38.0)**zzn)
c	write(8,*) i,fBB,fE,Fc,fSR(i),Fwca1
c	write(8,*) i,fBB,Fwpy,Fwsy,Fwgy,Fwcy,Fwpa,Fwsa,Fwga,Fwca,Fmp,
c     x    Fms,Fmg,Fmc
	 Fbp = (1/alphas)*((DLSOC(i)-dlsy)*Fwsy+(DLSOC(i)-dlsa)*Fwsa 
     x  +(DLSOC(i)-dlpy)*Fwpy+(DLSOC(i)-dlpa)*Fwpa +(DLSOC(i)-dlsa)*Fms 
     x  +(DLSOC(i)-dlpa)*Fmp)
        Fbg=(1/alphac(i))*((DLCOC(i)-dlcy)*Fwcy+(DLCOC(i)-dlca)*Fwca 
     x	+(DLCOC(i)-dlgy)*Fwgy+(DLCOC(i)-dlga)*Fwga +(DLCOC(i)-dlca)*Fmc 
     x  +(DLCOC(i)-dlga)*Fmg)
        Fbs=Fwpy+Fwpa+Fwsy+Fwsa +Fms +Fmp-Fbp
        Fbc=Fwgy+Fwga+Fwcy+Fwca +Fmc +Fmg-Fbg

        if(i.lt.571) O = O +(Fbg+(15./8.)*Fbp)*Dt -(Fwgy+Fwga+Fmg)*Dt 
     x   -(15./8.)*(Fwpy+Fwpa+Fmp)*Dt

        Spy=Spy+(Fbp-Fwpy-Fyop)*Dt
        Ssy =Ssy +(Fbs-Fwsy-Fyos)*Dt
        Gy= Gy+ (Fbg-Fwgy-Fyog)*Dt
        Cy=Cy + (Fbc-Fwcy-Fyoc)*Dt
c  CT is set near the start of the program
        Ca = CT-Gy - Ga - Cy -2.0
        dlpy = dlpy +((DLSOC(i)-dlpy-alphas)*Fbp/Spy)*Dt
        dlpa = dlpa + (Fyop*(dlpy-dlpa)/Spa)*Dt
        dlsy = dlsy+((DLSOC(i)-dlsy)*Fbs/Ssy)*Dt
        dlsa = dlsa + (Fyos*(dlsy-dlsa)/Ssa)*Dt
        dlgy=dlgy+((DLCOC(i)-dlgy-alphac(i))*Fbg/Gy)*Dt
        dlga =dlga+(Fyog*(dlgy-dlga)/Ga)*Dt
        dlcy=dlcy+((DLCOC(i)-dlcy)*Fbc/Cy)*Dt
        dlca=dlca+(Fyoc*(dlcy-dlca)/Ca)*Dt
        FSi0=6.67
        fAD(i)=fA(i)*fD(i)
        fB=(Fbc-(Fwcy+Fwca))/((fAD(i)**0.65)*fE*fR(i)*FSi0)
c	print *,fB,Fbc,Fwcy,Fwca,fAD(i),fE,fR(i),FSi0
c        print *,'GYM,LIFE,ACT,FERT'
        n=0
        if(i.gt.381) then 
          RCO2old=RCO2*2.
c	  if((i/10)*10.eq.i) print *,'in loop ',i
          do while(abs(RCO2/RCO2old-1.0).gt.0.001)
	    RCO2old=RCO2
            Fbbs= (RCO2**(0.5+ACT*GCM))*((1.0+RT*GCM*alog(RCO2)
     x        -RT*Ws*fac+RT*GEOG(i))**0.65)
     x        *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            W=((0.5+ACT*GCM)*RCO2**(-0.5+ACT*GCM))*((1+RT*GCM*alog(RCO2)
     x       -RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x       *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            V=(RCO2**(0.5+ACT*GCM))*0.65*((1.0+RT*GCM*alog(RCO2)
     x        -RT*Ws*(fac)+RT*GEOG(i))**(-0.35))*(RT*GCM/RCO2)
     x         *exp(-ACT*Ws*fac)*exp(ACT*GEOG(i))
            fPBS = W + V
            if(RCO2.gt.((Fbbs-fB)/fPBS)) then
c damp the iteration to avoid overshoot
              RCO2=RCO2-0.9*((Fbbs-fB)/fPBS)
	    else
c restrict the iteration to avoid nonpositive value
              RCO2=RCO2*0.2
	    endif
	  end do
  777 format(i5,8g13.5)
        elseif(i.le.381.and.i.gt.351) then
          RCO2old=RCO2*2.
	  n=0
          do while(abs(RCO2/RCO2old-1.0).gt.0.001)
            RCO2old=RCO2
            oldfBBS=(RCO2**(0.5+ACT*GCM))*((1+RT*GCM*alog(RCO2)
     x        -RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x        *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            oldW=(0.5+ACT*GCM)*RCO2**(-0.5+ACT*GCM)
     x          *((1.0+RT*GCM*alog(RCO2)-RT*Ws*(fac)
     x          +RT*GEOG(i))**0.65)
     x          *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            oldV=RCO2**(0.5+ACT*GCM)*0.65*((1.0+RT*GCM*alog(RCO2)
     x       -RT*Ws*(fac)+RT*GEOG(i))**(-0.35))*(RT*GCM/RCO2)
     x       *exp(-ACT*Ws*fac)*exp(ACT*GEOG(i))
            oldfPBS = oldW + oldV
	    oldX=0.
            ewfBBS=((2.0**FERT)*RCO2**(FERT+ACT*GCM))
     x       *((1+RCO2)**(-FERT))
     x   *((1.0+RT*GCM*alog(RCO2)-RT*Ws*fac+RT*GEOG(i))**0.65)
     x       *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
           ewW=(2.0**FERT)*(FERT+ACT*GCM)*(RCO2**(FERT+ACT*GCM-1.0))
     x 	    *((1.0+RCO2)**(-FERT))*((1.0+RT*GCM*alog(RCO2)
     x 	    -RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x 	    *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            ewV=(-FERT*(1.0+RCO2)**(-(1.0+FERT)))
     x	  *((2**FERT)*RCO2**(FERT+ACT*GCM))*((1.0+RT*GCM*alog(RCO2)
     x	  -RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x	  *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            ewX=0.65*((1.0+RT*GCM*alog(RCO2)-RT*Ws*(fac) 
     x     +RT*GEOG(i))**(-0.35))*(RT*GCM/RCO2)
     x	 *(2**FERT*RCO2**(FERT+ACT*GCM))*((1+RCO2)**(-FERT))
     x	 *exp(-ACT*Ws*fac)*exp(ACT*GEOG(i))
            ewfPBS=ewW+ewV+ewX
	    fBBS=((i-351)/30.)*oldfBBS + ((381-i)/30.)*ewfBBS
	    fPBS=((i-351)/30.)*oldfPBS + ((381-i)/30.)*ewfPBS
            if(RCO2.gt.((fBBS-fB)/fPBS)) then
c damp the iteration to avoid overshoot
              RCO2=RCO2-0.9*((fBBS-fB)/fPBS)
	    else
c restrict the iteration to avoid nonpositive value
              RCO2=RCO2*0.2
	    endif
            n=n+1
          end do
        else	
          RCO2old=RCO2*2.
c	  if((i/10)*10.eq.i) print *,'in loop ',i
          do while(abs(RCO2/RCO2old-1.0).gt.0.001)
	    RCO2old=RCO2
            Fbbs=((2**FERT)*RCO2**(FERT+ACT*GCM))*((1+RCO2)**(-FERT))
     x	         *((1+RT*GCM*alog(RCO2)-RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x	         *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            W=(2**FERT)*(FERT+ACT*GCM)*(RCO2**(FERT+ACT*GCM-1))
     x	       *((1+RCO2)**(-FERT))*((1+RT*GCM*alog(RCO2)
     x	       -RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x	       *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
       V=(-FERT*(1+RCO2)**(-(1.+FERT)))*((2**FERT)*RCO2**(FERT+ACT*GCM))
     x	    *((1+RT*GCM*alog(RCO2)-RT*Ws*(fac)+RT*GEOG(i))**0.65)
     x	    *exp(-ACT*Ws*(fac))*exp(ACT*GEOG(i))
            X=0.65*((1+RT*GCM*alog(RCO2)-RT*Ws*(fac) 
     x         +RT*GEOG(i))**(-0.35))*(RT*GCM/RCO2)
     x	       *(2**FERT*RCO2**(FERT+ACT*GCM))*((1+RCO2)**(-FERT))
     x	       *exp(-ACT*Ws*fac)*exp(ACT*GEOG(i))
            fPBS=W+V+X
            if(RCO2.gt.((Fbbs-fB)/fPBS)) then
c damp the iteration to avoid overshoot
              RCO2=RCO2-0.9*((Fbbs-fB)/fPBS)
	    else
c restrict the iteration to avoid nonpositive value
              RCO2=RCO2*0.2
	    endif
          end do 
        endif
        tau=15 + 6*alog(RCO2)-12.8*fac+GEOG(i)
        O2 =100*( O/(O+143.))
        ppm=250.*RCO2
c  the indexing here is time equals (i-1)Ma, save the ppm values at 10-Ma intervals, K=1 --> 0Ma
	if(((i-1)/10)*10.eq.i-1) then
	  Ma=-i+1
c          print *, Ma, ppm
c	  write(8,*) -i,ppm
          k=1+i/10
          pp(k1,k2,k3,k4,k)=ppm
	endif
        gcsppm(i)=ppm
	aage(i)=i-1
      end do
c  compute variance reduction -- compare variance of residual to variance of data
c  normalized -- sum(pred-data)^2/sum(data)^2
      sumsq=0.
      sumdatsq=0.
      do i=1,ndat
        diff=alog(gcsppm(iage(i))/pppm(i))
	sumsq=sumsq+(diff)**2
        diff=alog(pppm(i)/250.)
	sumdatsq=sumdatsq+(diff)**2
      end do
      varred=sumsq/sumdatsq
      pp(k1,k2,k3,k4,59)=varred
      pp(k1,k2,k3,k4,60)=varred
      varredmin=amin1(varredmin,pp(k1,k2,k3,k4,60))
	end do
	end do
	end do
c        print *,'GYM,LIFE,ACT,FERT'
        print *,GYM,LIFE,ACT,FERT
	end do
c      close(8)
c  sift through the runs to find acceptable data fits  (chiquare <= NDAT; bias < 0.3)
      do k=1,2
      do j=1,4
      do i=1,10
        jkk(i,j,k)=0
      end do
      end do
      end do
        do k1=1,10
        do k2=1,10
        do k3=1,10
        do k4=1,10
          if(pp(k1,k2,k3,k4,59).le.varredmin*1.5) then
		  jkk(k1,1,1)=jkk(k1,1,1)+1
		  jkk(k2,2,1)=jkk(k2,2,1)+1
		  jkk(k3,3,1)=jkk(k3,3,1)+1
		  jkk(k4,4,1)=jkk(k4,4,1)+1
	  endif
	  if(pp(k1,k2,k3,k4,60).le.varredmin*1.2) then
		  jkk(k1,1,2)=jkk(k1,1,2)+1
		  jkk(k2,2,2)=jkk(k2,2,2)+1
		  jkk(k3,3,2)=jkk(k3,3,2)+1
		  jkk(k4,4,2)=jkk(k4,4,2)+1
	  endif
        end do
        end do
        end do
        end do
	print *,'deltaT=',deltaT
	print *,'minimum chi-squared misfit',varredmin
	print *,'parameters for runs that fit with varred 50% from min'
	print *,'ACT'
	print *,(jkk(i,1,1),i=1,10)	  
	print *,'FERT'
	print *,(jkk(i,2,1),i=1,10)	  
	print *,'LIFE'
	print *,(jkk(i,3,1),i=1,10)	  
	print *,'GYM'
	print *,(jkk(i,4,1),i=1,10)	  
	print *,'   '
	print *,'parameters for runs that fit with varred 20% from min'
	print *,'ACT'
	print *,(jkk(i,1,2),i=1,10)	  
	print *,'FERT'
	print *,(jkk(i,2,2),i=1,10)	  
	print *,'LIFE'
	print *,(jkk(i,3,2),i=1,10)	  
	print *,'GYM'
	print *,(jkk(i,4,2),i=1,10)
c      	  pause
     
c  order the ppm values at each age to form median and percentile values
      print *,'enter ordering loop'
      print *,'we order the values, and compute medians and percentiles'
      print *,'also use pdf weights on parameters, downweight extremes'
c      do i=1,5
c        wt(i)=0.01+i*0.03
c	wt(11-i)=wt(i)
c      end do   
      do i=1,5
        wt(i)=0.02+(i-1)*0.04
	wt(11-i)=wt(i)
      end do   
      print *,(wt(i),i=1,10)   
      do k=59,59
        do k1=1,10
	  kk1=1000*(k1-1)
        do k2=1,10
	  kk2=kk1+100*(k2-1)
        do k3=1,10
	  kk3=kk2+10*(k3-1)
        do k4=1,10
	  kk4=kk3+k4
	  dum(kk4)=pp(k1,k2,k3,k4,k)
	  weight(kk4)=wt(k1)*wt(k2)*wt(k3)*wt(k4)
	end do
	end do
	end do
	end do
	do i=1,9999
	  test=dum(i)
	  testwt=weight(i)
	  kmn=i
	  do j=i+1,10000
	    if(dum(kmn).gt.dum(j)) kmn=j
	  end do
	  dum(i)=dum(kmn)
	  dum(kmn)=test
	  weight(i)=weight(kmn)
	  weight(kmn)=testwt
	end do
c  integrate the weights into a cumulative probability curve
        do i=2,10000
	  weight(i)=weight(i)+weight(i-1)
	end do
        print *,iiit,weight(10000)
	print *,deltaT
	p01(iiit)=dum(100)
	p10(iiit)=dum(1000)
	p33(iiit)=dum(3333)
	p50(iiit)=dum(5000)
	p66(iiit)=dum(6666)
	p90(iiit)=dum(9000)
	p99(iiit)=dum(9900)
	ama(iiit)=deltaT
	S0(1)=1.
	S0(2)=1.
	S0(3)=1.
c	if((k/5)*5.eq.k) then
c	  call plotit(S0,dum,dm,10000,' ',' ',' ',1,0,0.,0,21)
c	  call plotit(S0,weight,dm,10000,' ',' ',' ',1,0,0.,0,22)
c	endif
c  lets compute something different, the fraction of runs that achieve better than N% variance residual
        i=1
	if(dum(1).gt.0.15) then
	  q01(iiit)=0.0
	else 
	  do while(dum(i).lt.0.15)
	    i=i+1
	  end do
	  q01(iiit)=float(i)/10000.
	endif
	iii=i
	if(dum(iii).gt.0.2) then
	  q10(iiit)=0.0
	else 
	  do while(dum(i).lt.0.20)
	    i=i+1
	  end do
	  q10(iiit)=float(i)/10000.
	endif
	iii=i
	if(dum(iii).gt.0.25) then
	  q33(iiit)=0.0
	else 
	  do while(dum(i).lt.0.25)
	    i=i+1
	  end do
	  q33(iiit)=float(i)/10000.
	endif
	iii=i
	if(dum(iii).gt.0.3) then
	  q50(iiit)=0.0
	else 
	  do while(dum(i).lt.0.3)
	    i=i+1
	  end do
	  q50(iiit)=float(i)/10000.
	endif
      end do
      end do
      q66(1)=q01(1)
      q90(1)=q10(1)
      q99(1)=q33(1)
      do i=2,ndel
        q66(i)=q66(i-1)+q01(i)
        q90(i)=q90(i-1)+q10(i)
        q99(i)=q99(i-1)+q33(i)
      end do
      do i=1,ndel
        q66(i)=q66(i)/q66(ndel)
        q90(i)=q90(i)/q90(ndel)
        q99(i)=q99(i)/q99(ndel)
      end do
      open(8,file='out1.dat',form='formatted')
      write(8,106) (ama(i),p01(i),p10(i),p33(i),p50(i),p66(i),p90(i),
     x  p99(i),i=1,ndel)
      close(8)
      open(8,file='out2.dat',form='formatted')
      write(8,106) (ama(i),q01(i),q10(i),q33(i),q50(i),q66(i),q90(i),
     x  q99(i),i=1,ndel)
      close(8)
      write(title,105) deltaT
  105 format('misfit residual for CO\\sub{2} excess ',f3.1)
  106 format(8f12.4)
c      call plotit_axes(0.,0.,0.,1.0)
c      call plotit(ama,p01,dm,ndel,' ',' ',' ',2,0,0.03,0,0)
c      call plotit(ama,p99,dm,ndel,' ',' ',' ',2,0,0.03,0,0)
c      call plotit(ama,p10,dm,ndel,' ',' ',' ',2,0,0.06,0,0)
c      call plotit(ama,p90,dm,ndel,' ',' ',' ',2,0,0.06,0,0)
c      call plotit(ama,p33,dm,ndel,' ',' ',' ',2,0,0.1,0,0)
c      call plotit(ama,p66,dm,ndel,' ',' ',' ',2,0,0.1,0,0)
c      call plotit(ama,p50,dm,ndel,title,'DeltaT',
c     x 'fractional variance',2,0,0.0,0,1)
c      call plotit_axes(0.,0.,0.,0.0)
c      call plotit(ama,p01,dm,ndel,' ',' ',' ',2,2,0.03,0,0)
c      call plotit(ama,p99,dm,ndel,' ',' ',' ',2,2,0.03,0,0)
c      call plotit(ama,p10,dm,ndel,' ',' ',' ',2,2,0.06,0,0)
c      call plotit(ama,p90,dm,ndel,' ',' ',' ',2,2,0.06,0,0)
c      call plotit(ama,p33,dm,ndel,' ',' ',' ',2,2,0.1,0,0)
c      call plotit(ama,p66,dm,ndel,' ',' ',' ',2,2,0.1,0,0)
c      call plotit(ama,p50,dm,ndel,title,'DeltaT',
c     x 'fractional variance',2,2,0.0,0,1)
c      call plotit_axes(0.,0.,0.,1.0)
c      call plotit(ama,q01,dm,ndel,' ',' ',' ',2,0,0.03,0,0)
c      call plotit(ama,q99,dm,ndel,' ',' ',' ',2,0,0.03,0,0)
c      call plotit(ama,q10,dm,ndel,' ',' ',' ',2,0,0.06,0,0)
c      call plotit(ama,q90,dm,ndel,' ',' ',' ',2,0,0.06,0,0)
c      call plotit(ama,q33,dm,ndel,' ',' ',' ',2,0,0.1,0,0)
c      call plotit(ama,q66,dm,ndel,' ',' ',' ',2,0,0.1,0,0)
c      call plotit(ama,q50,dm,ndel,'weighted '//title,'DeltaT',
c     x 'fractional variance',2,0,0.0,0,1)
c      call plotit_axes(0.,0.,0.,0.0)
c      call plotit(ama,q01,dm,ndel,' ',' ',' ',2,2,0.03,0,0)
c      call plotit(ama,q99,dm,ndel,' ',' ',' ',2,2,0.03,0,0)
c      call plotit(ama,q10,dm,ndel,' ',' ',' ',2,2,0.06,0,0)
c      call plotit(ama,q90,dm,ndel,' ',' ',' ',2,2,0.06,0,0)
c      call plotit(ama,q33,dm,ndel,' ',' ',' ',2,2,0.1,0,0)
c      call plotit(ama,q66,dm,ndel,' ',' ',' ',2,2,0.1,0,0)
c      call plotit(ama,q50,dm,ndel,'weighted '//title,'DeltaT',
c     x 'fractional variance',2,2,0.0,0,1)
      stop
      end
