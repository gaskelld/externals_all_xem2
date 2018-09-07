C=======================================================================
                                                                        
      SUBROUTINE QUASIY8(E0,EP,TH,SIGMA)                                
                                                                        
C Calculates quasielastic cross section 
c using Donnelly/Sick super scaling model
c but Paris w.f. for deuteron
c for IA=1, does simple resolution smearing only
c changed to call F1F2QE 7/08 pyb

      IMPLICIT NONE     
      REAL E0,EP,TH,SIGMA
      REAL  CHBAR,ALPHA,PM,PI                                            
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828)                                       
      PARAMETER (PI    = 3.1415927)                                     
      REAL FF/1./  !phony initialization to please compiler: SER 4/15/93
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM 
      COMMON     /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      logical smrelasp, usearen
      real smrdEp,dep, dW
      common/smrproton/ smrelasp,smrdEp,usearen
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL PAULI_SUP1,PAULI_SUP2,FFF
      REAL*8 GEP,GEN,GMP,GMN
      real*8 qv,a,b,c,amf2,ams,disc,y,dydep,en,emp,denominator
      REAL THR,SINSQ,COSTH,QSQ,TAU,W1,W2,CMOTT,CSMOTT,RECOIL,DSIGE,FY,SD
      real*8 alfa,c1,c2,ARENHOVEL_SIG
      real*8 kappa,lam,lamp,taup,squigglef,psi,psip,nuL,nut
      real*8 kf,es,GM2bar,GE2bar,W1bar,W2bar,Delta,GL,GT
      common/testing/prttst,usegrd
      integer i,model,in_range
      logical prttst,usegrd
      real*8 z8,a8,q28,w28,f18,f28

! 100*prob(Pz) from Paris wave function in 10 MeV bins from -1 to 1 GeV
! integeral is 1. Integral is 100.
      integer izz
      real*8 depdpz,pz,qvp
! using just w*w for d-state (tail too big)
      real*8 fydr(200)/
     > 0.00000,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00000/
!*** NO D-sate (tail too small)
       real*8 fydn(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00030,
     > 0.00033,0.00038,0.00041,0.00047,0.00050,0.00057,0.00061,0.00068,
     > 0.00073,0.00081,0.00087,0.00095,0.00102,0.00111,0.00119,0.00128,
     > 0.00137,0.00147,0.00156,0.00167,0.00176,0.00189,0.00198,0.00210,
     > 0.00220,0.00232,0.00242,0.00254,0.00264,0.00276,0.00284,0.00296,
     > 0.00303,0.00313,0.00319,0.00327,0.00331,0.00338,0.00340,0.00344,
     > 0.00344,0.00347,0.00345,0.00348,0.00346,0.00349,0.00350,0.00358,
     > 0.00366,0.00386,0.00412,0.00454,0.00510,0.00594,0.00704,0.00859,
     > 0.01060,0.01334,0.01688,0.02156,0.02763,0.03544,0.04561,0.05878,
     > 0.07561,0.09741,0.12537,0.16167,0.20811,0.26892,0.34765,0.45096,
     > 0.58609,0.76539,1.00428,1.32177,1.75207,2.33167,3.11006,4.13194,
     > 5.42556,6.92172,8.36256,9.28786,9.28786,8.36256,6.92172,5.42556,
     > 4.13194,3.11006,2.33167,1.75207,1.32177,1.00428,0.76539,0.58609,
     > 0.45096,0.34765,0.26892,0.20811,0.16167,0.12537,0.09741,0.07561,
     > 0.05878,0.04561,0.03544,0.02763,0.02156,0.01688,0.01334,0.01060,
     > 0.00859,0.00704,0.00594,0.00510,0.00454,0.00412,0.00386,0.00366,
     > 0.00358,0.00350,0.00349,0.00346,0.00348,0.00345,0.00347,0.00344,
     > 0.00344,0.00340,0.00338,0.00331,0.00327,0.00319,0.00313,0.00303,
     > 0.00296,0.00284,0.00276,0.00264,0.00254,0.00242,0.00232,0.00220,
     > 0.00210,0.00198,0.00189,0.00176,0.00167,0.00156,0.00147,0.00137,
     > 0.00128,0.00119,0.00111,0.00102,0.00095,0.00087,0.00081,0.00073,
     > 0.00068,0.00061,0.00057,0.00050,0.00047,0.00041,0.00038,0.00033,
     > 0.00030,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/
! using 1.5 * (1-cos**2)**2 for dtate (better than either of above)
       real*8 fyd(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00031,
     > 0.00033,0.00038,0.00042,0.00048,0.00052,0.00058,0.00063,0.00070,
     > 0.00076,0.00084,0.00090,0.00099,0.00107,0.00117,0.00125,0.00136,
     > 0.00145,0.00158,0.00168,0.00181,0.00192,0.00207,0.00219,0.00235,
     > 0.00249,0.00265,0.00280,0.00297,0.00313,0.00332,0.00347,0.00368,
     > 0.00385,0.00406,0.00424,0.00447,0.00467,0.00491,0.00513,0.00541,
     > 0.00565,0.00597,0.00626,0.00665,0.00703,0.00751,0.00803,0.00868,
     > 0.00938,0.01030,0.01135,0.01265,0.01421,0.01616,0.01848,0.02143,
     > 0.02496,0.02942,0.03489,0.04170,0.05015,0.06060,0.07374,0.09018,
     > 0.11064,0.13651,0.16893,0.21017,0.26205,0.32886,0.41420,0.52472,
     > 0.66768,0.85544,1.10347,1.43058,1.87107,2.46124,3.25006,4.28232,
     > 5.58540,7.08967,8.53679,9.46493,9.46493,8.53679,7.08967,5.58540,
     > 4.28232,3.25006,2.46124,1.87107,1.43058,1.10347,0.85544,0.66768,
     > 0.52472,0.41420,0.32886,0.26205,0.21017,0.16893,0.13651,0.11064,
     > 0.09018,0.07374,0.06060,0.05015,0.04170,0.03489,0.02942,0.02496,
     > 0.02143,0.01848,0.01616,0.01421,0.01265,0.01135,0.01030,0.00938,
     > 0.00868,0.00803,0.00751,0.00703,0.00665,0.00626,0.00597,0.00565,
     > 0.00541,0.00513,0.00491,0.00467,0.00447,0.00424,0.00406,0.00385,
     > 0.00368,0.00347,0.00332,0.00313,0.00297,0.00280,0.00265,0.00249,
     > 0.00235,0.00219,0.00207,0.00192,0.00181,0.00168,0.00158,0.00145,
     > 0.00136,0.00125,0.00117,0.00107,0.00099,0.00090,0.00084,0.00076,
     > 0.00070,0.00063,0.00058,0.00052,0.00048,0.00042,0.00038,0.00033,
     > 0.00031,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/
      logical doing_elas,usef1f2qe
      common/doelas/ doing_elas

      SIGMA = 0.                                                        
      IF (iA.LE.1 .and. (.not.smrelasp)) RETURN
      IF (iA.LE.1 .and. (.not.doing_elas)) RETURN

      if(IA.eq.2 .and. usearen .and. (.not.doing_elas)) then
cc        model = 1
        model = 2
        sigma = ARENHOVEL_SIG(E0,EP,TH,MODEL,IN_RANGE)
        return
      endif

c xxx changed to use new routine unless arenhovel specified for D2
      usef1f2qe = .true.
      if(usef1f2qe .and. (.not.doing_elas)) then
       THR   = TH*PI/180.                                                
       SINSQ = SIN(THR/2.)**2                                            
       COSTH = COS(THR)                                                  
       QSQ   = 4.*E0*EP*SIN(THR/2.)**2                                   
       z8 = iz
       a8 = ia
       q28 = qsq
       w28 = pm**2 + 2. * pm * (e0 - ep) - qsq
 
c       CALL  F1F2QE07(z8, a8, q28, w28, F18, F28) 
c       CALL  F1F2QE09(z8, a8, q28, w28, F18, F28) 
       CALL  F1F2QE16(z8, a8, q28, w28, F18, F28) 
  
       W2 = F28 / (e0 - ep)
       W1 = F18 / PM
       CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
       CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
       sigma  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT
       if(abs(e0-5.157).lt.0.001.and.abs(ep-3.55).lt.0.01.and.
     >   abs(th-23.).lt.0.01)
     >  write(6,'(''dbgqe'',8f8.3)') e0,ep,th,q28,w28,F18,F28,sigma
       return
      endif

      DYDEP = 0.                                                        
      AMS   = avgM-PM                                                   
      THR   = TH*PI/180.                                                
      SINSQ = SIN(THR/2.)**2                                            
      COSTH = COS(THR)                                                  
      QSQ   = 4.*E0*EP*SIN(THR/2.)**2                                   
      TAU   = QSQ/4.0/PM**2                                             
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                          
c      if(prttst) write(8,'(1x,i3,11f7.3)') 
c     >  ig,e0,ep,th,qsq,gep,gmp,gmn,gen

                                                                        
! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     >  PAULI_SUP1,PAULI_SUP2)
      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)                           
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU)  
c     if(prttst) write(8,'(10f7.3)') pauli_sup2,gep**2,gen**2,w1,tau,w2
                                                                        
      CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
      CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
      RECOIL = PM/(PM+E0*(1.-COSTH))                                    
      DSIGE  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT*RECOIL                 

! new for proton. Use dw = e * dep / ep = dep / recoil
      if(doing_elas) then
        if(avgA.gt.1.5) then
         call NUC_FORM_FACTOR(QSQ,W1,W2,FFF)
         RECOIL = PM * avgA / (PM * avgA + E0 * (1.-COSTH)) 
         DSIGE  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT*RECOIL                 
cc         write(6,'(''dbg elas 1'',5f10.3)') QSQ,W1,W2,FFF,recoil
        endif
        dep = EP - E0 * recoil
        sigma = 0.
cc        write(6,'(''dbg elas 2'',5f10.3)') e0,ep,dep,smrdEp
        if(abs(dep/smrdEp).lt.5.) then
          sigma = dsige / smrdep  / sqrt(3.1416) * 
     >      exp(-dep**2 / smrdEp**2)
          if(sigma.lt.0.) write(6,'(''dbg elas 2'',7f8.3,3e11.4)') 
     >     e0,ep,dep,smrdEp,recoil,w1,w2,csmott,dsige,sigma
        endif
        return
      endif

! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
      if(IA.eq.2) kf=0.085
      if(iA.eq.2) Es=0.0022
      if(IA.eq.3) kf=0.180
      if(iA.eq.3) Es=0.010 
c Use narrowe with for q.e.
cc      if(IA.eq.4) kf=0.200
      if(IA.eq.4) kf=0.180
      if(iA.eq.4) Es=0.015 
cc      if(iA.eq.4) Es=0.012 
      if(IA.gt.4) kf=0.220
      if(iA.gt.4) Es=0.015 
      if(IA.gt.7) kf=0.225
      if(iA.gt.7) Es=0.020 
      if(IA.gt.14) kf=0.225
      if(iA.gt.14) Es=0.020 
      if(IA.gt.16) kf=0.230
      if(iA.gt.16) Es=0.025 
      if(IA.gt.25) kf=0.236
      if(iA.gt.25) Es=0.018 
      if(IA.gt.38) kf=0.241
      if(iA.gt.38) Es=0.028 
      if(IA.gt.55) kf=0.241
      if(iA.gt.55) Es=0.023 
      if(IA.gt.60) kf=0.245
      if(iA.gt.60) Es=0.028 
      if(ep.ge.e0) then 
cc        write(6,'(1x,''error,e0,ep='',2f8.3)') e0,ep
        return
      endif
      QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)                          
      kappa = qv / 2. / pm
      lam = (e0 - ep) / 2. / pm
      if(abs(kappa**2 - lam**2 - tau).gt.0.01) then
        write(6,'(1x,''error,tau='',3f8.3)') kappa**2, lam**2,tau
        return
      endif
      lamp = lam - Es / 2. / pm
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/pm)**2) -1.
      if(1.+lamp.le.0.) then
        write(6,'(1x,''error,lamp='',3f8.3)') lam,lamp
        return
      endif
      if(taup * (1. + taup).le.0.) then
        write(6,'(1x,''error,taup='',3f8.3)') kappa**2, lam**2,tau
        return
      endif
      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))
      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))
      nuL = (tau / kappa**2)**2
      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
      GM2bar = Pauli_sup1 * (iZ * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (iZ * GEP**2 + avgN * GEN**2) 
      W1bar = tau * GM2bar
      W2bar = (GE2bar + tau * GM2bar) / (1. + tau)
      Delta = squigglef * (1. - psi**2) * (
     >  sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >  (1. - psi**2) * tau / kappa**2)
      GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
      GT = (2. * tau * GM2bar + Delta * W2bar) /
     >  2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)

c added to prevent negative xsections:
      gt = max(0., gt)

! Correction factors that make better agreement with world data
! on He and C from Day and NE5. See cn15.f
! BUT these don't work for Q2<0.3, so take them out again.
! actually, much better off without them as far as he/C concerned
      if(IA.gt.2.and.ia.lt.8) then
cc        GL = GL * (1.577 - 0.485 / qv)                
cc        GT = GT * (0.875 + 0.494 / qv)                   
      endif
      if(ia.ge.8) then
cc        GL = GL * (1.875 - 0.286 / qv)                
cc        GT = GT * (1.022 + 0.235 / qv)                   
      endif

! from Maria Barbaro: see superfy.m1
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
      sigma = csmott * FY * (nuL * GL + nuT * GT) / kf 
      if(sigma.lt.0.) then
cc        write(6,'(''ERROR, sigma...='',10f8.1)')
cc     >    sigma,csmott,fy,nuL,GL,nuT,GT
        sigma = 0.
      endif
      if(prttst) write(8,'(/1x,12f7.3)') e0,ep,th,W1,W2,csmott/10000.,
     >  fy,nul,gl,nut,gt,sigma/100.

! Use PWIA and Paris W.F. for deuteron
      if(IA.eq.2) then
! add extra effecive width by makeing qv wider for resolution
        qvp = qv
        if(smrelasp) qvp = qv * sqrt(1. + 
     >    (1.0 * smrdEp / 2. / 0.055 / qv)**2)
! qvp never used: THIS IS NOT IMPLEMENTED!!!!!
        pz = (2.*pm*(e0-ep) - qsq) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        depdpz = ep * qv / pm / e0
! Facotr of 100 and 0.01 cance to give:
        sigma = dsige * fyd(izz) / depdpz  
        if(prttst) write(8,'(i4,4f7.3)') izz,pz,depdpz,
     >    fyd(izz),sigma/100.
      endif

      return

! Are we changing quasi-elastic cross section due to Asymmetry?         
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0)                                
c    >     .AND.(INDEX(TARGET,'_P').GT.0) )                             
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >  CSMOTT*RECOIL,DSIGE)                                            
                                                                        
      QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)                          
      AMF2  = PM**2                                                     
      A     = 4.*QV**2-4.*(E0-EP+avgM)**2                               
      B     = -4.*QV*(AMS**2-AMF2-QV**2+(E0-EP+avgM)**2)                
      C     = (E0-EP+avgM)**4+(AMS**2-AMF2-QV**2)**2                    
     >       -2.*(E0-EP+avgM)**2*(AMF2+QV**2+AMS**2)                    
      DISC  = B**2-4.*A*C                                               
      IF (A.EQ.0..OR.DISC.LT.0.) RETURN                                 
      Y     = (-B-SQRT(DISC))/2./A                                      
                                                                        
C JACOBIAN DY/DEP                                                       
                                                                        
      EN    = SQRT(AMF2+(QV+Y)**2)                                      
      EMP   = SQRT(AMS**2+Y**2)                                         
      DENOMINATOR = (Y/EMP+(QV+Y)/EN)                                   
      IF (DENOMINATOR.LE.1.E-30) RETURN                                 
      DYDEP = (1.+(QV+Y)/EN*(EP-E0*COSTH)/QV)/DENOMINATOR               
                                                                        
C GAUSSIAN SMEARING FUNCTION                                            
                                                                        
      FY = 0.                                                           
      IF (iA.LE.4) THEN                                                 
! *** 034 way too small for 4He: swith to 100 
c**          SD = 0.034                                         
           SD = 0.100
      ELSE                                                              
           SD = 0.119                                                   
      END IF                                                            
c special case for C
c***      if(IA.eq.12) SD=0.09

      fy=0.
      IF ((Y**2/2./SD**2).LT.40.)                                       
     >     FY = 1./SD/SQRT(2.*PI)*EXP(-Y**2/2./SD**2)                   

! Modified to use fit to FY from nucl-th/9812078 (C. delgi Atti...)
! Modified 11/05 to get better agreement with Naida Fromin's data
      if(IA.eq.2) then
        a = 6.1 - 0.50
        alfa = .045 
        c1 = 0.018 
        c2 = 0.25 
      endif
      if(IA.eq.3) then
        a = 7.1 - 1.00
        alfa = .083
        c1 = 0.041 * 1.00
        c2 = 0.33 * 1.20
      endif
      if(IA.eq.4) then
        a = 6.8 - 1.00
        alfa = .167
        c1 = 0.106 * 1.20
        c2 = 0.65 * 1.20
      endif
      if(IA.gt.4.and.IA.le.30) then
        a = 5.1
        alfa = .166
        c1 = 0.083
        c2 = 0.57 * 1.20
      endif
      if(IA.gt.30) then
        a = 4.6
        alfa = .138
        c1 = 0.058
        c2 = 0.62 * 1.30
      endif
!*** peb pyb changed exponent from -6 to -8 10/05
!*** to try to get better agreement with data
! changed to 8 on 11/15/05
      fy = c1 * exp(-1.0 * a**2 * y**2) / (alfa**2 + y**2) +
     >   c2 * exp(-8.0 * abs(y))

      SIGMA = DSIGE*FY*DYDEP                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                           
