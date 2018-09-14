C=======================================================================
                                                                        
      SUBROUTINE F1F2IN09(Z, A, QSQ, Wsq, F1, F2)                       
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
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 Z,A,qsq,wsq,w,f1c,x
      real*8 avgn,r,dr,nu,eps,kappa,sigres,flux,siglp,sigtp,F1pp,F1dp
      REAL*8 W1,W2,sigt,rc,w1pb,w2pb,F1,F2,sigl,F1d, F1p,qv
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des,fyuse
      real*8 xp,f1cor
      real a4, x4, fitemc, emcfac
      logical goodfit
      INTEGER ISM

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

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

cc     
       real*8 F1M
       logical DEBUG/.TRUE./

      IA = int(A)
      avgN = A - Z
      nu = (wsq - pm**2 + qsq) / 2. / pm
      qv = sqrt(nu**2 + qsq)

      if(Wsq.le.0.0) W = 0.0
      W  = sqrt(Wsq)
      x  = QSQ / (2.0 * pm * nu)
      if(Wsq.le.0.0) x = 0.0

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
! changed 4/09
        if(IA.eq.3) kf=0.115
        if(iA.eq.3) Es=0.001 
! changed 4/09
        if(IA.gt.3) kf=0.19
        if(iA.gt.3) Es=0.017
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
c changed 5/09
        if(iA.gt.7) Es=0.0165

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.230
c

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
! changed 5/09 
        if(iA.gt.55) Es=0.018 
 

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
          xp = qsq/(wsqp-pm*pm+qsq)

          IF(WSQP.GT. 1.159) THEN
            call CHRISTY507(Wsqp,Qsq,F1pp,Rc,sigtp,siglp)
            call resmodd(wsqp,qsq,xvald0,F1dp)
c            f1cor = 1.0+xvalc(26)*xp+xvalc(27)*xp**2+xvalc(28)*xp**3+
c     &          xvalc(29)*xp**4+xvalc(30)*xp**5

            f1cor = xvalc(26)*(1.0-xvalc(27)*xp*xp)**xvalc(28)*
     &               (1.0-xvalc(29)*exp(-1.0*xvalc(30)*xp))
           

c            f1cor = f1cor/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq))

c            write(6,*) qsq,xp,f1cor
            F1dp = F1dp*f1cor
            F1pp = F1pp*f1cor
            F1d = F1d + F1dp * Fyuse
            F1p = F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse   !!!  Not ready to use yet
c            call rescsp(Wsqp,Qsq,sigtp,siglp)
c            call rescsn(Wsqp,Qsq,sigtn,sigln)
c            F1d =  F1d+(F1n + F1pp) * Fyuse
c            F1p =  F1p + F1pp * Fyuse
            sigt = sigt + sigtp * Fyuse
            sigl = sigl + siglp * Fyuse
          ENDIF
        ENDDO
        Rc = 0.
        if(sigt .gt. 0.) Rc = sigl / sigt
        W1 = (2. * Z * F1d + (A - 2. * Z) * (2. * F1d - F1p)) / PM 


CC TEST BELOW if commented out CC


c        W1= W1*(1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)
c        write(6,*) P(13),P(14),P(15),P(16),P(17)
c        W1= W1*(1.0+xvalc(26)*x+xvalc(27)*x**2+xvalc(28)*x**3+
c     &          xvalc(29)*x**4+xvalc(30)*x**5)

cc        if(W .GT. 1.3) 

c         if(W .GT. 0.0) 
c     >       W1=W1*(1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ))**2

        CALL MEC2009( Z , A , qsq , wsq , F1M )

        W1 = W1 + F1M
c        if(Wsq .gt.0.0 ) Rc = Rc * ( 1.0 + P(6) + P(23)*A )
        if(Wsq .gt.0.0 ) Rc = Rc * (1.0 + xvalc(35)*x + 
     &          xvalc(36)*x*x/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq)))
c        Rc = Rc/(1.0+xvalc(31)*log(1.0+xvalc(32)*qsq))

        W2 = W1 * (1. + Rc) / (1. + nu**2 / qsq)

        DEBUG=.FALSE.
        IF( W1 .LE. 0.0 .AND. DEBUG ) THEN 
           write(*,*) 'test  = ', Z,A,W,QSQ,x,F1M,W1
           write(*,*) 'test1 = ', 
     >          (1.0+(P(20)*W+P(21)*W**2)/(1.0+P(22)*QSQ)),
     >        (1.0+P(13)*x+P(14)*x**2+P(15)*x**3+P(16)*x**4+P(17)*x**5)
        ENDIF

      ENDIF

      A4 = A
      x4 = qsq / 2. / pm / nu
      emcfac = fitemc(x4, a4, goodfit)

      F1 = pm * W1 
      F2 = nu * W2  
c      F1 = pm * W1 * emcfac 
c      F2 = nu * W2 * emcfac 



      RETURN                                                            
      END                                          

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MEC2009(z,a,q2,w2,f1)

! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,am/0.9383/,w,nu
      real*8 a1,b1,c1,t1,dw2,dw22,emc,w2min,damp
      real*4 x4,a4,fitemc
      logical goodfit
      integer i
      real*8 pb(20)/ 
     >     0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >     0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >    -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >     0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >     0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      real*8 p18

      real*8 x, f1corr 

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

      f1 = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - am**2 + q2) / 2. / am
      x  = q2 / (2.0 * am * nu )

      if(a.lt.2.5) return

      p18 = p(18)
! special case for 3He
      if(a.gt.2.5 .and. a.lt.3.5) p18 = 70
! special case for 4He
      if(a.gt.3.5 .and. a.lt.4.5) p18 = 170.
! new values for C, Al, Cu
      if(a.gt.4.5) p18 = 215.
      if(a.gt.20.) p18 = 235.
      if(a.gt.50.) p18 = 230.

c      write(6,*) "here2 ", xvalc(1),xvalc(2),xvalc(3),xvalc(5)

       
       f1corr = P(0)*exp(-((W-P(1))**2)/(P(2)))/ 
     >      ((1.0 + MAX( 0.3 , Q2 ) / P(3) ) ** P(4) )*nu**P(5)
     >      *( 1.0 + P18 * A ** ( 1.0 + P(19) * x ) )

       if(a.eq.12) then

        a1 = q2**2*xvalc(1)*exp(-1.0*q2/xvalc(2))/
     &                                (xvalc(3)+q2)**xvalc(4)

        b1 = xvalc(5)+xvalc(6)*q2

c        c1 = xvalc(7)/(1.+xvalc(8)*sqrt(q2))

        c1 = xvalc(33)+xvalc(34)*q2

c         c1 = 0.290   !!! Test  !!!

        t1 = (w2-b1)**2/c1**2/2.
c        dw2 = w2+q2/xvalc(10)-1.*am*am
c        dw22 = w2+q2/3.-1.*am*am
        dw2 = w2+q2*(1.-1./2.2)-1.0*am*am

        if(dw2.LT.0.0) dw2 = 0.0
c        if(dw22.LT.0.0) dw22 = 0.0
        f1corr = a1*(exp(-1.*t1)*sqrt(dw2))

        f1corr = f1corr+xvalc(38)*exp(-1.*((w2-1.25)/xvalc(39))**2.0)*
     &               q2*exp(-1.*xvalc(40)*q2)


        x4 = x
        a4 = a
        goodfit = .true.
c        emc = fitemc(x4,a4,goodfit)
c        f1corr = f1corr/am/emc
        f1corr = f1corr/am
       endif

       f1 = f1corr

       if(f1 .le.1.0E-9 ) f1=0.0
c       write(*,*) 'vahe1= ', A, W*W, Q2, f1corr

      return
      end
