CCC  Version 040207  -  Author:  M.E. Christy                       CCC
CCC  This routine returns proton photo-absorbtion cross sections     CCC
CCC  for either transverse or longitudinal photons in units of       CCC
CCC  ub/Sr/Gev.                                                      CCC
CCC  Fit form is empirical.  Interpret physics from it at your       CCC
CCC  own risk.           
C=======================================================================
                                                                        
      SUBROUTINE F1F2IN07(Z, A, QSQ, Wsq, F1, F2)                       
!--------------------------------------------------------------------
! Fit to inelastic cross sections for A(e,e')X
! valid for all W<3 GeV and all Q2<10 GeV2
! 
! Inputs: Z, A (real*8) are Z and A of nucleus 
!         (use Z=0., A=1. to get free neutron)
c         Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c         Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
! Version of 10/20/2006 P. Bosted
!--------------------------------------------------------------------
      implicit none
      real*8 Z,A,qsq,wsq,w,f1c,x
      real*8 avgn,r,dr,nu,eps,kappa,sigres,flux,siglp,sigtp,F1pp,F1dp
      REAL*8 W1,W2,sigt,rc,w1pb,w2pb,F1,F2,sigl,F1d, F1p,qv
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des,fyuse
      real a4, x4, fitemc, emcfac
      logical goodfit
      INTEGER ISM

! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/

! This is for exp(-xx**2/2.), from teste.f
       real*8 xxp(99)/
     > -3.000,-2.939,-2.878,-2.816,-2.755,-2.694,-2.633,-2.571,-2.510,
     > -2.449,-2.388,-2.327,-2.265,-2.204,-2.143,-2.082,-2.020,-1.959,
     > -1.898,-1.837,-1.776,-1.714,-1.653,-1.592,-1.531,-1.469,-1.408,
     > -1.347,-1.286,-1.224,-1.163,-1.102,-1.041,-0.980,-0.918,-0.857,
     > -0.796,-0.735,-0.673,-0.612,-0.551,-0.490,-0.429,-0.367,-0.306,
     > -0.245,-0.184,-0.122,-0.061, 0.000, 0.061, 0.122, 0.184, 0.245,
     >  0.306, 0.367, 0.429, 0.490, 0.551, 0.612, 0.673, 0.735, 0.796,
     >  0.857, 0.918, 0.980, 1.041, 1.102, 1.163, 1.224, 1.286, 1.347,
     >  1.408, 1.469, 1.531, 1.592, 1.653, 1.714, 1.776, 1.837, 1.898,
     >  1.959, 2.020, 2.082, 2.143, 2.204, 2.265, 2.327, 2.388, 2.449,
     >  2.510, 2.571, 2.633, 2.694, 2.755, 2.816, 2.878, 2.939, 3.000/
! these are 100x bigger for convenience
       real*8 fyp(99)/
     > 0.0272,0.0326,0.0390,0.0464,0.0551,0.0651,0.0766,0.0898,0.1049,
     > 0.1221,0.1416,0.1636,0.1883,0.2159,0.2466,0.2807,0.3182,0.3595,
     > 0.4045,0.4535,0.5066,0.5637,0.6249,0.6901,0.7593,0.8324,0.9090,
     > 0.9890,1.0720,1.1577,1.2454,1.3349,1.4254,1.5163,1.6070,1.6968,
     > 1.7849,1.8705,1.9529,2.0313,2.1049,2.1731,2.2350,2.2901,2.3379,
     > 2.3776,2.4090,2.4317,2.4454,2.4500,2.4454,2.4317,2.4090,2.3776,
     > 2.3379,2.2901,2.2350,2.1731,2.1049,2.0313,1.9529,1.8705,1.7849,
     > 1.6968,1.6070,1.5163,1.4254,1.3349,1.2454,1.1577,1.0720,0.9890,
     > 0.9090,0.8324,0.7593,0.6901,0.6249,0.5637,0.5066,0.4535,0.4045,
     > 0.3595,0.3182,0.2807,0.2466,0.2159,0.1883,0.1636,0.1416,0.1221,
     > 0.1049,0.0898,0.0766,0.0651,0.0551,0.0464,0.0390,0.0326,0.0272/

      integer iz,ia,i
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/

! deuteron fit parameters
       real*8 xvald0(50)/
     >  0.1964E+01, 0.1086E+01, 0.5313E-02, 0.1265E+01, 0.8000E+01,
     >  0.2979E+00, 0.1354E+00, 0.2200E+00, 0.8296E-01, 0.9578E-01,
     >  0.1094E+00, 0.3794E+00, 0.8122E+01, 0.5189E+01, 0.3290E+01,
     >  0.1870E+01, 0.6110E+01,-0.3464E+02, 0.9000E+03, 0.1717E+01,
     >  0.4335E-01, 0.1915E+03, 0.2232E+00, 0.2119E+01, 0.2088E+01,
     > -0.3029E+00, 0.2012E+00, 0.1104E-02, 0.2276E-01,-0.4562E+00,
     >  0.2397E+00, 0.1204E+01, 0.2321E-01, 0.5419E+03, 0.2247E+00,
     >  0.2168E+01, 0.2266E+03, 0.7649E-01, 0.1457E+01, 0.1318E+00,
     > -0.7534E+02, 0.1776E+00, 0.1636E+01, 0.1350E+00,-0.5596E-02,
     >  0.5883E-02, 0.1934E+01, 0.3800E+00, 0.3319E+01, 0.1446E+00/
                                                                        
      IA = int(A)
      avgN = A - Z
      nu = (wsq - pm**2 + qsq) / 2. / pm
      qv = sqrt(nu**2 + qsq)

! Cross section for proton or neutron
      W1 = 0.
      W2 = 0.
      IF(IA .lt. 2 .and. wsq.gt.1.155) THEN
        call CHRISTY507(Wsq,Qsq,F1p,Rc,sigt,sigl)
! If neutron, subtract proton from deuteron. Factor of two to
! convert from per nucleon to per deuteron
        if(Z .lt. 0.5) then
          call resmodd(wsq,qsq,xvald0,F1d)
          F1p = F1d * 2.0 - F1p
        endif
        W1 = F1p / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      ENDIF

! For deuteron
      if(IA .eq. 2) then
c get Fermi-smeared R from Erics proton fit
        call pind(Wsq, Qsq, F1c, Rc, sigt, sigl)
c get fit to F1 in deuteron, per nucleon
        call resd(qsq, wsq, xvald0, F1d)
c convert to W1 per deuteron
        W1 = F1d / PM * 2.0
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)
      endif

! For nuclei
      IF(IA.gt.2) then
        sigt = 0.
        sigl = 0.
        F1d = 0.
        F1p = 0.
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
       if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
! changed 9/24/07
        if(IA.eq.3) kf=0.100
        if(iA.eq.3) Es=0.010 
! changed 9/24/07
        if(IA.eq.4) kf=0.170
        if(iA.eq.4) Es=0.015 
        if(IA.gt.4) kf=0.165
        if(iA.gt.4) Es=0.015 
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
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
cc for test
        if(iA.gt.60) Es=0.018

! adjust pf to give right width based on kf
        pf = 0.5 * kf 
! assume this is 2 * pf * qv
        DW2DPF = 2. * qv
        dw2des = 2. * (nu + PM) 
! switched to using 99 bins!
cc        do ism = 1,15
cc          fyuse = fy(ism)
cc          WSQP = WSQ + XX(ISM) * PF * DW2DPF - es * dw2des
        do ism = 1,99
          fyuse = fyp(ism)/100.
          WSQP = WSQ + XXp(ISM) * PF * DW2DPF - es * dw2des
          IF(WSQP.GT. 1.159) THEN
            call CHRISTY507(Wsqp,Qsq,F1pp,Rc,sigtp,siglp)
            call resmodd(wsqp,qsq,xvald0,F1dp)
            F1d = F1d + F1dp * Fyuse
            F1p = F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse
          ENDIF
        ENDDO
        Rc = 0.
        if(sigt .gt. 0.) Rc = sigl / sigt
        W1 = (2. * Z * F1d + (A - 2. * Z) * (2. * F1d - F1p)) / PM 
        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq) 
      ENDIF

      A4 = A
      x4 = qsq / 2. / pm / nu
      emcfac = fitemc(x4, a4, goodfit)

! this version is per nucleon
      F1 = pm * W1 * emcfac  / A
      F2 = nu * W2 * emcfac  / A

! Add MEC correction

      RETURN                                                            
      END                                                               

CCC  Version 051407  -  Author:  M.E. Christy                               CCC
CCC  This fit version includes data from E00-116 (see thesis of S. Malace)  CCC 
CCC  as well preliminary data from E00-002.                                 CCC  
CCC  Subroutine to get Transvese and Longitudinal eP cross sections         CCC 
CCC  from fits cross sections over a range of epsilon.  The subroutine      CCC
CCC  resmod.f is required.  Units are in ub/Sr/Gev.                         CCC


      SUBROUTINE christy0507(W2,Q2,F1,R)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end

CCC  Version 040207  -  Author:  M.E. Christy                       CCC
CCC  This routine returns proton photo-absorbtion cross sections     CCC
CCC  for either transverse or longitudinal photons in units of       CCC
CCC  ub/Sr/Gev.                                                      CCC
CCC  Fit form is empirical.  Interpret physics from it at your       CCC
CCC  own risk.                                                       CCC


      SUBROUTINE christy507(W2,Q2,F1,R,sigt,sigl)

      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),xval(100)
      real*8 mp,mp2,pi,alpha,xb,sigT,sigL,F1,FL,F2,R
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036


      data xval / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /


      do i=1,50
        xval1(i) = xval(i)
        xvalL(i) = xval(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
 
 
      xb = q2/(w2+q2-mp2)

      call resmod507(1,w2,q2,xval1,sigT)
      call resmod507(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = sigL/sigT
    
      end
      SUBROUTINE RESMOD507(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      INTEGER i,j,l,num,sf


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
           
    
      dip = 1./(1.+q2/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.91)**rescoef(i,4)
 
        else

          height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &                             *exp(-1.*rescoef(i,3)*q2)

        endif

 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
        height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1) then

        do i=1,2  

          h_nr(i) = nr_coef(i,1)/     
     &       (q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        sig_nr = sig_nr*xpr(1)


      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       


 1000  format(8f12.5)

      RETURN 
      END 

      subroutine pind(W2,Q2,F1,R,sigt,sigl)
! Calculate proton with Fermi smearing of a deuteron 
      implicit none
      real*8 q2,w2,F1,R,sigt,sigl,am/0.9383/,nu,qv,F1p,Rp,sigtp,siglp
      real*8 amd/1.8756/,w2p,pz
      integer ism
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
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
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1=0.
      R=0.
      sigt=0.
      sigl=0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
         w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fyd(ism) / 10.
          sigl = sigl + siglp * fyd(ism) / 10.
          F1   = F1   + F1p   * fyd(ism) / 10.
        endif
       enddo
      else
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
c Need avp2f term to get right behavior x>1! 
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call CHRISTY507(W2p,Q2,F1p,Rp,sigtp,siglp)
          sigt = sigt + sigtp * fydf(ism) / 100.
          sigl = sigl + siglp * fydf(ism) / 100.
          F1   = F1   + F1p   * fydf(ism) / 100.
        endif
       enddo
      endif

      if(sigt.ne.0.) R = sigl / sigt
      return
      end
      
      subroutine resd(q2,w2,xval,F1)
! Calculate dueteron F1 by Fermi smearing of proton plus neutron 
! Add MEC term (not smeared)
      implicit none
      real*8 q2,w2,xval(50),F1,am/0.9383/,nu,qv,dw2dpf,w2p,sigp,f1sv
      real*8 sigtst,amd/1.8756/, pz, f1m,f2m
      integer ism,i
       real*8 fyd(20)/ 0.4965, 0.4988, 0.4958, 0.5008, 0.5027, 0.5041,
     >  0.5029, 0.5034, 0.4993, 0.5147, 0.5140, 0.4975, 0.5007, 0.4992,
     >  0.4994, 0.4977, 0.5023, 0.4964, 0.4966, 0.4767/
       real*8 avpz(20)/-0.1820,-0.0829,-0.0590,-0.0448,-0.0345,-0.0264,
     > -0.0195,-0.0135,-0.0079,-0.0025, 0.0029, 0.0083, 0.0139, 0.0199,
     >  0.0268, 0.0349, 0.0453, 0.0598, 0.0844, 0.1853/
       real*8 avp2(20)/ 0.0938, 0.0219, 0.0137, 0.0101, 0.0081, 0.0068,
     >  0.0060, 0.0054, 0.0051, 0.0049, 0.0050, 0.0051, 0.0055, 0.0060,
     >  0.0069, 0.0081, 0.0102, 0.0140, 0.0225, 0.0964/
c Look up tables for deuteron in fine bins for sub threshold
       real*8 fydf(200)/
     > 0.00001,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
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
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00001/
       real*8 avp2f(200)/
     >     1.0,0.98974,0.96975,0.96768,0.94782,0.94450,0.92494,0.92047,
     > 0.90090,0.89563,0.87644,0.87018,0.85145,0.84434,0.82593,0.81841,
     > 0.80021,0.79212,0.77444,0.76553,0.74866,0.73945,0.72264,0.71343,
     > 0.69703,0.68740,0.67149,0.66182,0.64631,0.63630,0.62125,0.61154,
     > 0.59671,0.58686,0.57241,0.56283,0.54866,0.53889,0.52528,0.51581,
     > 0.50236,0.49291,0.47997,0.47063,0.45803,0.44867,0.43665,0.42744,
     > 0.41554,0.40656,0.39511,0.38589,0.37488,0.36611,0.35516,0.34647,
     > 0.33571,0.32704,0.31656,0.30783,0.29741,0.28870,0.27820,0.26945,
     > 0.25898,0.25010,0.23945,0.23023,0.21943,0.20999,0.19891,0.18911,
     > 0.17795,0.16793,0.15669,0.14667,0.13553,0.12569,0.11504,0.10550,
     > 0.09557,0.08674,0.07774,0.06974,0.06184,0.05484,0.04802,0.04203,
     > 0.03629,0.03129,0.02654,0.02247,0.01867,0.01545,0.01251,0.01015,
     > 0.00810,0.00664,0.00541,0.00512,0.00512,0.00541,0.00664,0.00810,
     > 0.01015,0.01251,0.01545,0.01867,0.02247,0.02654,0.03129,0.03629,
     > 0.04203,0.04802,0.05484,0.06184,0.06974,0.07774,0.08674,0.09557,
     > 0.10550,0.11504,0.12569,0.13553,0.14667,0.15669,0.16793,0.17795,
     > 0.18911,0.19891,0.20999,0.21943,0.23023,0.23945,0.25010,0.25898,
     > 0.26945,0.27820,0.28870,0.29741,0.30783,0.31656,0.32704,0.33571,
     > 0.34647,0.35516,0.36611,0.37488,0.38589,0.39511,0.40656,0.41554,
     > 0.42744,0.43665,0.44867,0.45803,0.47063,0.47997,0.49291,0.50236,
     > 0.51581,0.52528,0.53889,0.54866,0.56283,0.57241,0.58686,0.59671,
     > 0.61154,0.62125,0.63630,0.64631,0.66182,0.67149,0.68740,0.69703,
     > 0.71343,0.72264,0.73945,0.74866,0.76553,0.77444,0.79212,0.80021,
     > 0.81841,0.82593,0.84434,0.85145,0.87018,0.87644,0.89563,0.90090,
     > 0.92047,0.92494,0.94450,0.94782,0.96768,0.96975,0.98974,1.0/
      logical usemec/.true./
c      common/mecc/usemec

      nu = (w2 - am**2 + q2) / 2. / am
      qv = sqrt(nu**2 + q2)
      F1 = 0.
! Do fast 20 bins if abvoe threshold
      if(w2.gt.1.30) then
       do ism = 1,20
        w2p = (amd + nu - sqrt(am**2 + avp2(ism)))**2 - 
     >    qv**2 + 2. * qv * avpz(ism) - avp2(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fyd(ism) / 10.
        endif
       enddo
      else
       f1 = 0.
       do ism = 1,200
        pz = -1. + 0.01 * (ism-0.5)
! Need avp2f term to get right behavior x>1!
        w2p = (amd + nu - sqrt(am**2 + avp2f(ism)))**2 - 
     >    qv**2 + 2. * qv * pz - avp2f(ism)
        if(w2p.gt.1.155) then
          call resmodd(w2p,q2,xval,sigp)
          F1 = F1 + sigp * fydf(ism) / 100. 
        endif
       enddo
      endif
! add MEC term
! took out again, then put back in again
      call mec(1.D0,2.D0,q2,w2,f1m,f2m,xval)
      if(usemec) f1 = f1 + f1m
      return
      end

CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC
CCC *** This is modified to give fit to deuteron for sigt
CCC *** Modified to be for sigt only (former sf=1)
CCC *** Modified to take out lowq2 business.
CCC *** params 1-12 are now hard-wired and used for n/p ratio
ccc *** instead
ccc added code to pre-calculate w-dependent parameters
ccc differences from Eric 507
c a) definition of A**2 (factor 2 M W)
c b) powers of L in width for 2-pion
c c) 0.71 -> 0.91 in dipole formula
             
! changed to use F1 instead of sigt
      SUBROUTINE RESMODD(w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(5000,7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,xpr,alpha,pi,F1
      INTEGER i,j,l,lmax,num,iw
      real*8 noverp,fnp_nmc,x,a,b,sig_mec,brp(7,3)
      real*8 xval0(12)/
c new 5/07 values. Values 1 and 7 will be overridden below.
     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00/
      real*8 xvalold(50),w2sv,q2sv,sigrsv(7),md,w2p,wp,wdifp,xprp,nu
      logical first/.true./
      common/tst2/sigrsv,sig_nr,sig_mec
      real br2(7),br3(7)
      save

      sig = 0.
      if(w2.lt.1.07327**2 .or. w2.gt.25 .or. 
     >  q2.lt.0.0 .or. q2.gt.11.0) then
        write(6,'(1x,''error, q2 or w2 out of range'',2f8.3)') w2,q2
        return
      endif

c do this if fitting masses or widths, else set first true in above
      first=.false.
      if(xvalold(1).ne.xval(1) .or.
     >   xvalold(7).ne.xval(7)) first=.true.
      xvalold(1)=xval(1)
      xvalold(7)=xval(7)

      if(first) then
       mp = 0.9382727
       mpi = 0.135
       mpi2 = mpi*mpi
       meta = 0.547
       mp2 = mp*mp
       pi = 3.141593
       alpha = 1./137.036

! branching ratios
       br(1,1) = 1.0     
       br(2,1) = 0.5
       br(3,1) = 0.65
       br(4,1) = 0.65
       br(5,1) = 0.4
       br(6,1) = 0.65
       br(7,1) = 0.6

! angular momenta
       ang(1) = 1.       !!!  P33(1232)
       ang(2) = 0.       !!!  S11(1535)
       ang(3) = 2.       !!!  D13(1520)
       ang(4) = 3.       !!!  F15(1680)
       ang(5) = 0.       !!!  S15(1650)
       ang(6) = 1.       !!!  P11(1440) roper   
       ang(7) = 3.       !!!  ? 4th resonance region

! x0 parameter
       do i=1,7
c 2006
c         x0(i) = 0.165
c 2007
         x0(i) = 0.215
       enddo
c 2006
c      x0(4) = 0.6
c 2007
       x0(1) = 0.1446

! out branching ratio
       do i=1,7
         br(i,2) = 1.-br(i,1)
       enddo
    
! remember w2
       w2sv = w2

! uses xvals of 1-12, 47, and 48
! move masses, wdiths into local variables
! pyb changed to be fixed
       num = 0
       do i=1,6              
         num = num + 1
         mass(i) = xval0(i)
       enddo
       do i=1,6             
         num = num + 1
         intwidth(i) = xval0(num)
       enddo
! changed to allow delta width, mass to vary
! taken out again since xval(1) used in MEC
c       mass(1) = xval(1)
c       intwidth(1) = xval(7)
c 2006
c       mass(7) = xval(47)
c       intwidth(7) = xval(48)
c 2007
       mass(7) = 1.9341
       intwidth(7) = 0.380

! precalculate w-dependent quantites in 0.1 MeV bins
       do iw=1073,5000
        w = 0.001 * (iw+0.5)
        w2 = w**2
        wdif = w - (mp + mpi)
        wr = wdif/w

! Calculate kinematics needed for threshold Relativistic B-W 
        k = (w2 - mp2) / 2. / mp
        kcm = (w2 - mp2) / 2. / w
        epicm = (W2 + mpi**2 -mp2 ) / 2. / w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 ) / 2. / w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 ) / 2. / w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))
        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))
! Calculate partial widths
          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
          if(i.ne.2) then
c            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
c     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)
c make same as resmod507
            pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     >       **(ang(i)+2.)
     >       * W / mass(i)
          else
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif
          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)
          pgam(i) = intwidth(i)*pgam(i)
          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(iw,i) = width(i) * pgam(i) / ((W2 - mass(i)**2.)**2. 
     &            + (mass(i)*width(i))**2.) *
     >            kr(i) / k * kcmr(i) / kcm / intwidth(i)
        enddo ! loop on i
c        write(55,'(i5,f7.2,7f10.4)') iw,w,(sigr(iw,i),i=1,7)
       enddo ! loop on iw

       w2 = w2sv
       first = .false.
       do i=1,50
         xvalold(i) = xval(i)
       enddo
! write table for article
c       open(unit=9,file='nhd.tbl1')
       do i=1,7
         br2(i)=br(i,2)
         br3(i)=0.
         if(i.eq.2) br2(i)=0.
         if(i.eq.2) br3(i)=br(i,2)
         if(i.le.6) then
c          write(9,199) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),(xval(12 + (i-1)*4 + j),j=1,4)
c 199      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & ',f7.1,' & ',f6.3,' & ',f6.3,' \\\\')
         else
c          write(9,198) i,mass(i),intwidth(i),
c     >     int(ang(i)),X0(i),
c     >     br(i,1),br2(i),br3(i),xval(49)
c 198      format(i1,' & ',f6.3,' & ',f6.3,' & ',i1,' & ',f5.3,
c     >      ' & ',f4.2,' & ',f4.2,' & ',f4.2,' & ',f6.3,
c     >      ' & 0.0 & 0.0 & 4.0 \\\\')
         endif
       enddo
c       close(unit=9)
c       open(unit=9,file='nhd.tbl2')
c       write(9,197) (xval(i),i=2,6)
       do i=1,2
c         write(9,197) (xval(36 + (i-1)*4 + j),j=1,4),xval(44+i)
c 197     format(f7.1,' & ',f7.4,' & ',f7.4,' & ',f7.4,' & ',
c     >     f7.4,' \\\\')
       enddo
c       write(9,'(''xval50='',f10.4)') xval(50)
       close(unit=9)
      endif ! if first
      
! get parameters into local variables
      num = 12
! resonance height coefficients. xvals of 13-36
      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo
      enddo
!  Non-Res coefficients xvals of 37-44
      do i=1,2               
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo

! Begin resonance Q^2 dependence calculations   CCC
! uses xvals 49
      do i=1,6
        height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2) * q2 / (1. + rescoef(i,3) * q2))/
ccc     &          (1. + q2/0.71)**rescoef(i,4)
c make same as resmod507
     &          (1. + q2/0.91)**rescoef(i,4)
      enddo
ccc      dip = 1./(1. + q2 / 0.71)**2  
c make same as resmod507
      dip = 1./(1. + q2 / 0.91)**2  
      dip2 = dip**2
ccc      height(7) = xval(49)*dip2 
c make same as resmod507
      height(7) = xval(49)*dip 
      iw = int(1000.*sqrt(w2))
      sig_res = 0.
      do i=1,7
ccc        sigrsv(i) =  height(i) * sigr(iw,i)
c make same as resmod507 by squaring height
        sigrsv(i) =  height(i)**2 * sigr(iw,i)
        sig_res = sig_res + sigrsv(i) 
      enddo

c make same as resmod507
      sig_res = sig_res * sqrt(w2)

! Begin non-resonant part uses xvals 45, 46, 50
! Depends on both W2 and Q2 so can't easily precalculate
      sig_nr = 0.
c      xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
! to make same as resmod507
      xpr = 1.+(w2-(mp+mpi)**2)/(q2+0.05)
      xpr = 1./xpr
      w = sqrt(w2)
      wdif = w - (mp + mpi)
      do i=1,2  
        sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &   /(q2+nr_coef(i,2))**
     &   (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
      enddo
c new section for sqrt(W)
c make same as resmod507 by turning this off
ccc        sig_nr = sig_nr +(xval(2)*(wdif)**0.5)
ccc     &   /(q2 + xval(3))**
ccc     &   (xval(4) + xval(5) * q2 + xval(6) * q2**2)
        
      sig_nr = sig_nr * xpr
     
! Add third term to try to describe MEC, now using Wdiff in 
! deuteron rather than proton
! ** Taken out 10/17/06 
c     md = 2.*mp
c     nu = (q2 + w2 - mp2) / 2. / mp
c     w2p = md**2 + 2. * md * nu - q2
c     Wp = sqrt(w2p)
c     wdifp = wp - md
c     sig_mec = 0.
c     if(wdifp .gt. 0.) then
c       xprp = 1. + (w2p - (md)**2) / (q2 + xval(50))
c       xprp = 1. / xprp
c       sig_mec = (xval(1) + xval(2)*wdifp**(1/2.) +
c    >     xval(3)*wdifp) /
c    &    (q2 + xval(4))**(xval(5) + xval(6) * q2) *
c    >    xprp
c      endif
c     sig = sig_res + sig_nr + sig_mec

      sig = sig_res + sig_nr 
c      write(6,'(1x,i4,2f7.2,4e10.3)') iw,q2,w2,height(1),
c     >  sigr(iw,1),sig_res,sig_nr

! changed to use F1 instead of sigt
      F1 = sig * (w2-mp2)/8./pi/pi/alpha/0.3894e3
      sig = F1

      RETURN 
      END 

      subroutine mec(z,a,q2,w2,f1,f2,xval)
! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,f2,xval(50),am/0.9383/,w,nu

      w = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am

! changed to use max(0.3,q2)
      f1 = xval(1) * exp(-(w - xval(2))**2/xval(3)) /
     >   (1. + max(0.3,q2) / xval(4))**xval(5) * nu ** xval(6)
      
      f2 = 0.
      if(q2.gt.0.) f2 = nu * (f1/am) / (1. + nu**2 / q2) 
      return
      end
