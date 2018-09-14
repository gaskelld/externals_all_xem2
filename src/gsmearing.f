C=======================================================================
                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, F1, F2)
                       
!--------------------------------------------------------------------

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,f1mec,f2mec
      real*8 nu,x,mp,mp2,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at
      real*8 xxp(100),fytot,fytot2,norm

      real a4, x4
      real*8 fitemct,emcfac,emcfacL,xfr,xfl
      logical goodfit
      INTEGER ISM,drn,wfn
      external off_mKP_fit

      real*8 xvalc(40) /
     & 0.42590E+00,0.84169E+01,0.34460E+00,0.68082E+01,0.81354E+00,
     & 0.12912E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97981E+00,0.96986E+00,0.10395E+01,0.99471E+00,0.98433E+00,
     & 0.10000E+01,0.97860E+00,0.10060E+01,0.98912E+00,0.99476E+00,
     & 0.99330E+00,0.99866E+00,0.10000E+01,0.10105E+01,0.10077E+01,
     & 0.10780E+01,0.48919E+00,0.12220E+01,0.13706E+00,0.41167E+02,
     & 0.60174E+01,0.12541E-01,0.17868E+00,0.54704E-01,-.22587E+00,
     & 0.18227E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+0 /

      drn = 5
      xfr = 0.95
      xfl = 1.E-3
      wfn = 2

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)
      a4 = A                                                    
      x4 = x

      if(A.GE.12.) then
        deltae = 0.015  !!! energy shift !!!
        kf = 0.228      !!! fermi momentum  !!!
      elseif(A.GE.27) then
        deltae = 0.017
        kf = 0.238
      elseif(A.GE.56) then
        deltae = 0.023
        kf = 0.241
      endif

      norm = 20.471
      norm = norm*2.0
c      norm = 1.0

      f1p = 0.0D0
      f1n = 0.0D0
      f2p = 0.0D0
      f2n = 0.0D0
      fLp = 0.0D0
      fLn = 0.0D0


c        sigt = 0.0D0
c        sigL = 0.0D0
      fytot = 0.0D0
      fytot2 = 0.0D0


! adjust pf to give right width based on kf
      pf = 0.5 * kf 
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + mp) 

      DO ism = 1,99

CCC   
        xxp(ism) = -3.0+6.0*float(ism-1)/98.0

        fyuse = 1.0D0/norm*exp(-0.5*xxp(ism)*xxp(ism))    !!! Gaussian !!!

CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des

CCC

        fytot = fytot+fyuse

c           write(6,2000) w2,q2,ism,xxp(ism),fyuse, wsqp, fytot


        IF(WSQP.GT. 1.159) THEN
          xt = q2/(q2+Wsqp-mp2)
          delta = off_mKP_fit(xt,wfn,drn)
          if(xt.GT.xfr) delta =  off_mKP_fit(xfr,wfn,drn)
          if(xt.LT.xfl) delta =  off_mKP_fit(xfl,wfn,drn)
          if(q2.LT.0.01) delta = 0.0

          offshell = 1.0D0/(1.0D0-delta)  

          offshell = 1.0D0      !!!  test

          x4 = xt

CCC   Next is medium modification factor2  CCC


          emcfac = xvalc(26)*(1.0-xvalc(27)*xt*xt)**xvalc(28)*
     &               (1.0-xvalc(29)*exp(-1.0*xvalc(30)*xt))

          emcfacL = (1.0 + xvalc(35)*xt + 
     &          xvalc(36)*xt*xt)/(1.0+xvalc(31)*log(1.0+xvalc(32)*q2)) 



          call sf(wsqp,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)
    
          f1pp = f1pp*emcfac*offshell
          f1nn = f1nn*emcfac*offshell
          fLpp = fLpp*emcfac*emcfacL*offshell
          fLnn = fLnn*emcfac*emcfacL*offshell
          f2pp = (2.*xt*f1pp+fLpp)/(1.+4.*xt*xt*mp2/q2)
          f2nn = (2.*xt*f1nn+fLnn)/(1.+4.*xt*xt*mp2/q2)

          F1p = F1p + F1pp * Fyuse
          F1n = F1n + F1nn * Fyuse
          F2p = F2p + F2pp * Fyuse
          F2n = F2n + F2nn * Fyuse      
          FLp = FLp + FLpp * Fyuse
          FLn = FLn + FLnn * Fyuse

         ENDIF

      ENDDO

      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

      call MEC2016(Z,A,w2,q2,f1mec)
      F1 = F1 + f1mec
      f2mec = 2.*x*f1mec/(1.+4.*x*x*mp2/q2)
      F2 = F2 + f2mec

c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MEC2016(z,a,w2,q2,f1corr)

! fit to low q2 dip region: purefly empirical
! assume contribution is pure transverse
      implicit none
      real*8 z,a,q2,w2,f1,mp/0.938272/,mp2,w,nu
      real*8 a1,b1,c1,t1,dw2,dw22,emc,w2min,damp
      integer i
      real*8 x, f1corr

      real*8 xvalc(40) / 
     & 0.42590E+00,0.84169E+01,0.34460E+00,0.68082E+01,0.81354E+00,
     & 0.12912E+00,0.15100E+01,0.16871E+01,0.33630E+00,-.22727E+01,
     & 0.97981E+00,0.96986E+00,0.10395E+01,0.99471E+00,0.98433E+00,
     & 0.10000E+01,0.97860E+00,0.10060E+01,0.98912E+00,0.99476E+00,
     & 0.99330E+00,0.99866E+00,0.10000E+01,0.10105E+01,0.10077E+01,
     & 0.10780E+01,0.48919E+00,0.12220E+01,0.13706E+00,0.41167E+02,
     & 0.60174E+01,0.12541E-01,0.17868E+00,0.54704E-01,-.22587E+00,
     & 0.18227E+01,0.10000E+01,0.00000E+00,0.10000E+01,0.00000E+0 /

      mp2 = mp*mp
      f1corr = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu )

      if(a.lt.2.5) return

      a1 = A*q2**2*xvalc(1)*exp(-1.0*q2/xvalc(2))/
     &                                (xvalc(3)+q2)**xvalc(4)

      b1 = xvalc(5)+xvalc(6)*q2

c        c1 = xvalc(7)/(1.+xvalc(8)*sqrt(q2))

      c1 = xvalc(33)+xvalc(34)*q2

c         c1 = 0.290   !!! Test  !!!


      t1 = (w2-b1)**2/c1**2/2.

      dw2 = w2+q2*(1.-1./2.2)-1.0*mp2

      if(dw2.LT.0.0) dw2 = 0.0
      f1corr = a1*(exp(-1.*t1)*sqrt(dw2))

      f1corr = f1corr+xvalc(38)*exp(-1.*((w2-1.25)/xvalc(39))**2.0)*
     &               q2*exp(-1.*xvalc(40)*q2)


c       f1corr = f1corr

       if(f1corr.LE.1.0E-9 ) f1corr=0.0


      return
      end
C=======================================================================
C ***********************************************************************
        FUNCTION off_mKP_fit (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]), constrained to vanish at x=0
C
C  Defined such that F2d = F2d(conv) + del^off F2d 
C       with off_mKP_fit = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Fit by Eric Christy (6/14/12).
C  Fortran code by W. Melnitchouk (6/17/12).
C
C  wfn: 1 (AV18)
C	2 (CD-Bonn)
C	3 (WJC-1)
C	4 (WJC-2)
C
C  dRN: 0 (0.0%)	[% change in nucleon radius in the deuteron]
C	1 (0.3%)
C	2 (0.6%)
C	3 (0.9%)
C	4 (1.2%)
C	5 (1.5%)
C	6 (1.8%)
C	7 (2.1%)
C	8 (2.4%)
C	9 (2.7%)
C      10 (3.0%)
C
C ***********************************************************************
	IMPLICIT NONE
	REAL*8	off_mKP_fit, x
	INTEGER	wfn, dRN
	REAL*8  p(0:9)

	off_mKP_fit = 0.D0
	IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
	IF (dRN.LT.0 .OR. dRN.GT.10) RETURN
! .......................................................................
	IF (wfn.EQ.1) THEN		! AV18
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02628D0
	    p(1) = 1.20029D0
	    p(2) = 7.49503D0
	    p(3) = 2.01901D0
	    p(4) = 0.00789D0
	    p(5) = 0.46739D0
	    p(6) = 0.73242D0
	    p(7) = 0.00328D0
	    p(8) = 0.87228D0
	    p(9) = 0.06400D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = 0.03638D0
	    p(1) = 0.38307D0
	    p(2) = 8.01156D0
	    p(3) = 2.30992D0
	    p(4) = 0.09027D0
	    p(5) = 0.69521D0
	    p(6) = 0.75973D0
	    p(7) = -0.05098D0
	    p(8) = 1.18963D0
	    p(9) = -0.19192D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02260D0
	    p(1) = 1.45377D0
	    p(2) = 0.50628D0
	    p(3) = 13.92200D0
	    p(4) = 0.03558D0
	    p(5) = 0.75147D0
	    p(6) = 0.86335D0
	    p(7) = -0.01383D0
	    p(8) = 1.04749D0
	    p(9) = 0.42099D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06410D0
	    p(1) = 1.18883D0
	    p(2) = 6.96799D0
	    p(3) = 8.87113D0
	    p(4) = 0.02603D0
	    p(5) = 0.70504D0
	    p(6) = 1.44139D0
	    p(7) = 0.00004D0
	    p(8) = -1.14305D0
	    p(9) = 0.73785D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06237D0
	    p(1) = 2.03192D0
	    p(2) = 4.01755D0
	    p(3) = 6.83741D0
	    p(4) = 0.04701D0
	    p(5) = -0.00457D0
	    p(6) = 1.30967D0
	    p(7) = -0.00996D0
	    p(8) = -0.42418D0
	    p(9) = 0.27524D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.06759D0
	    p(1) = 1.95103D0
	    p(2) = 3.54215D0
	    p(3) = 11.77533D0
	    p(4) = 0.09269D0
	    p(5) = 0.56534D0
	    p(6) = 0.98398D0
	    p(7) = -0.03031D0
	    p(8) = 3.26913D0
	    p(9) = -0.45923D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.07007D0
	    p(1) = 2.30938D0
	    p(2) = 4.94226D0
	    p(3) = 8.95701D0
	    p(4) = 0.06933D0
	    p(5) = 0.07145D0
	    p(6) = 1.94887D0
	    p(7) = -0.01210D0
	    p(8) = 5.92311D0
	    p(9) = 0.14312D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11965D0
	    p(1) = 2.06149D0
	    p(2) = 5.38881D0
	    p(3) = 12.08265D0
	    p(4) = 0.19668D0
	    p(5) = 0.61820D0
	    p(6) = 0.80489D0
	    p(7) = -0.08735D0
	    p(8) = 3.74802D0
	    p(9) = -0.70773D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.14735D0
	    p(1) = 2.27109D0
	    p(2) = 8.23092D0
	    p(3) = 7.31581D0
	    p(4) = 0.11953D0
	    p(5) = 0.67459D0
	    p(6) = 1.59118D0
	    p(7) = -0.02700D0
	    p(8) = 4.52840D0
	    p(9) = -1.77765D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.27194D0
	    p(1) = 2.01340D0
	    p(2) = 10.71380D0
	    p(3) = 8.84886D0
	    p(4) = 0.09345D0
	    p(5) = 0.49802D0
	    p(6) = 1.28523D0
	    p(7) = -0.00474D0
	    p(8) = 0.58703D0
	    p(9) = 0.88354D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.69848D0
	    p(1) = 1.48173D0
	    p(2) = 17.44991D0
	    p(3) = 12.73730D0
	    p(4) = 0.13118D0
	    p(5) = 0.34598D0
	    p(6) = 1.65884D0
	    p(7) = -0.02215D0
	    p(8) = 1.21306D0
	    p(9) = 0.96399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.2) THEN		! CD-Bonn
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02820D0
	    p(1) = 0.85879D0
	    p(2) = 9.48856D0
	    p(3) = 2.18885D0
	    p(4) = 0.00070D0
	    p(5) = -5.61817D0
	    p(6) = 14.80512D0
	    p(7) = 0.00348D0
	    p(8) = -1.30292D0
	    p(9) = -0.73075D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.02996D0
	    p(1) = 0.35717D0
	    p(2) = 6.53843D0
	    p(3) = 3.88389D0
	    p(4) = 0.00758D0
	    p(5) = -16.50399D0
	    p(6) = 77.60083D0
	    p(7) = 0.00320D0
	    p(8) = 0.42334D0
	    p(9) = 0.28545D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03261D0
	    p(1) = 0.91185D0
	    p(2) = 8.49348D0
	    p(3) = 10.19681D0
	    p(4) = 0.01598D0
	    p(5) = 0.83748D0
	    p(6) = 1.55960D0
	    p(7) = 0.00085D0
	    p(8) = -0.63447D0
	    p(9) = 0.65632D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.03034D0
	    p(1) = 1.58677D0
	    p(2) = 3.21753D0
	    p(3) = 11.66572D0
	    p(4) = 0.04999D0
	    p(5) = 0.56688D0
	    p(6) = 0.94941D0
	    p(7) = -0.01453D0
	    p(8) = -0.89157D0
	    p(9) = 0.27160D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.04831D0
	    p(1) = 1.75241D0
	    p(2) = 4.74662D0
	    p(3) = 8.29052D0
	    p(4) = 0.04730D0
	    p(5) = 0.33550D0
	    p(6) = 1.18790D0
	    p(7) = -0.00678D0
	    p(8) = -0.42800D0
	    p(9) = 0.36573D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.09019D0
	    p(1) = 1.22091D0
	    p(2) = 1.30114D0
	    p(3) = 17.58701D0
	    p(4) = 0.08312D0
	    p(5) = 0.66902D0
	    p(6) = 0.60767D0
	    p(7) = -0.02035D0
	    p(8) = 0.95978D0
	    p(9) = 1.11322D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.17060D0
	    p(1) = 1.42115D0
	    p(2) = 7.24672D0
	    p(3) = 5.80680D0
	    p(4) = 0.09200D0
	    p(5) = 0.43367D0
	    p(6) = 1.56378D0
	    p(7) = -0.02338D0
	    p(8) = 0.44968D0
	    p(9) = 0.29678D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11026D0
	    p(1) = 1.85213D0
	    p(2) = 6.74413D0
	    p(3) = 7.74362D0
	    p(4) = 0.08467D0
	    p(5) = 0.24708D0
	    p(6) = 1.12274D0
	    p(7) = -0.01505D0
	    p(8) = 0.44209D0
	    p(9) = 0.36126D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.15291D0
	    p(1) = 1.83333D0
	    p(2) = 7.76495D0
	    p(3) = 7.04783D0
	    p(4) = 0.09206D0
	    p(5) = 0.08655D0
	    p(6) = 1.27460D0
	    p(7) = -0.01659D0
	    p(8) = 0.45536D0
	    p(9) = 0.29407D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.24143D0
	    p(1) = 1.50401D0
	    p(2) = 9.33393D0
	    p(3) = 11.62779D0
	    p(4) = 0.09454D0
	    p(5) = 0.36361D0
	    p(6) = 0.82058D0
	    p(7) = -0.00802D0
	    p(8) = 0.34851D0
	    p(9) = 0.50844D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.22196D0
	    p(1) = 1.87228D0
	    p(2) = 10.18898D0
	    p(3) = 9.21038D0
	    p(4) = 0.11850D0
	    p(5) = 0.34360D0
	    p(6) = 1.28278D0
	    p(7) = -0.01754D0
	    p(8) = 0.54540D0
	    p(9) = 0.53457D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.3) THEN		! WJC-1
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.02322D0
	    p(1) = 0.11213D0
	    p(2) = 3.71079D0
	    p(3) = 5.51496D0
	    p(4) = 0.00877D0
	    p(5) = 0.84639D0
	    p(6) = 0.66227D0
	    p(7) = -0.00621D0
	    p(8) = -0.39896D0
	    p(9) = 0.32012D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.00058D0
	    p(1) = 2.33827D0
	    p(2) = 2.35664D0
	    p(3) = 36.75823D0
	    p(4) = -0.00752D0
	    p(5) = 0.05286D0
	    p(6) = 1.27262D0
	    p(7) = 0.01269D0
	    p(8) = 1.72720D0
	    p(9) = 0.20652D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03373D0
	    p(1) = 0.93858D0
	    p(2) = 0.15704D0
	    p(3) = 10.71630D0
	    p(4) = -0.00235D0
	    p(5) = -0.11937D0
	    p(6) = 0.74925D0
	    p(7) = 0.00452D0
	    p(8) = 2.96830D0
	    p(9) = -2.89070D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.08982D0
	    p(1) = 0.73060D0
	    p(2) = -0.16543D0
	    p(3) = 12.37035D0
	    p(4) = 0.04407D0
	    p(5) = 0.47361D0
	    p(6) = 0.74570D0
	    p(7) = -0.00933D0
	    p(8) = 0.53186D0
	    p(9) = 0.26943D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.11990D0
	    p(1) = 1.19824D0
	    p(2) = 3.06386D0
	    p(3) = 8.55017D0
	    p(4) = 0.05815D0
	    p(5) = 0.06123D0
	    p(6) = 1.45024D0
	    p(7) = -0.01414D0
	    p(8) = 0.48172D0
	    p(9) = 0.25171D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.15292D0
	    p(1) = 1.01991D0
	    p(2) = 1.20661D0
	    p(3) = 13.31860D0
	    p(4) = 0.02571D0
	    p(5) = -1.56438D0
	    p(6) = 2.69042D0
	    p(7) = -0.00000D0
	    p(8) = 0.29759D0
	    p(9) = 0.97967D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.35935D0
	    p(1) = 0.44637D0
	    p(2) = -0.25510D0
	    p(3) = 16.70057D0
	    p(4) = 0.10634D0
	    p(5) = 0.61659D0
	    p(6) = 0.58524D0
	    p(7) = -0.03335D0
	    p(8) = 0.93904D0
	    p(9) = 0.89819D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.97384D0
	    p(1) = -0.24934D0
	    p(2) = -0.61349D0
	    p(3) = 18.43254D0
	    p(4) = 0.18772D0
	    p(5) = 0.49599D0
	    p(6) = 0.61366D0
	    p(7) = -0.08116D0
	    p(8) = 0.87175D0
	    p(9) = 0.24026D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.21641D0
	    p(1) = 1.74710D0
	    p(2) = 5.19387D0
	    p(3) = 10.61285D0
	    p(4) = 0.06655D0
	    p(5) = 0.01300D0
	    p(6) = 0.94503D0
	    p(7) = -0.00642D0
	    p(8) = 0.48859D0
	    p(9) = 0.16331D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.32283D0
	    p(1) = 1.71708D0
	    p(2) = 7.51556D0
	    p(3) = 9.68202D0
	    p(4) = 0.09871D0
	    p(5) = 0.18788D0
	    p(6) = 0.80490D0
	    p(7) = -0.01673D0
	    p(8) = 0.48879D0
	    p(9) = 0.21016D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.32064D0
	    p(1) = 2.07514D0
	    p(2) = 9.34847D0
	    p(3) = 8.17225D0
	    p(4) = 0.10772D0
	    p(5) = 0.50272D0
	    p(6) = 1.30663D0
	    p(7) = -0.01215D0
	    p(8) = 0.59432D0
	    p(9) = 0.65399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.4) THEN		! WJC-2
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.03490D0
	    p(1) = 0.78902D0
	    p(2) = -0.25256D0
	    p(3) = 7.98679D0
	    p(4) = 0.00913D0
	    p(5) = 0.74835D0
	    p(6) = 0.60145D0
	    p(7) = -0.00464D0
	    p(8) = 0.41358D0
	    p(9) = 0.22524D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.01119D0
	    p(1) = 0.50514D0
	    p(2) = 19.35710D0
	    p(3) = 3.32395D0
	    p(4) = 0.00670D0
	    p(5) = 1.38279D0
	    p(6) = 1.24216D0
	    p(7) = 0.00049D0
	    p(8) = 0.38623D0
	    p(9) = 0.23497D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02653D0
	    p(1) = 1.27315D0
	    p(2) = -0.53410D0
	    p(3) = 14.08029D0
	    p(4) = 0.01474D0
	    p(5) = 1.82129D0
	    p(6) = 1.99455D0
	    p(7) = -0.00090D0
	    p(8) = 3.96583D0
	    p(9) = 4.61316D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06301D0
	    p(1) = 1.10373D0
	    p(2) = -0.26356D0
	    p(3) = 15.04038D0
	    p(4) = 0.02428D0
	    p(5) = -0.15349D0
	    p(6) = 3.03168D0
	    p(7) = 0.00127D0
	    p(8) = -0.73818D0
	    p(9) = 0.07474D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06150D0
	    p(1) = 2.15792D0
	    p(2) = 2.18241D0
	    p(3) = 9.84713D0
	    p(4) = 0.03608D0
	    p(5) = -0.13604D0
	    p(6) = 1.12241D0
	    p(7) = -0.00695D0
	    p(8) = -0.35646D0
	    p(9) = 0.31793D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.07179D0
	    p(1) = 1.97917D0
	    p(2) = 3.47662D0
	    p(3) = 10.00224D0
	    p(4) = 0.04587D0
	    p(5) = 0.06416D0
	    p(6) = 1.10677D0
	    p(7) = -0.00391D0
	    p(8) = -0.42677D0
	    p(9) = 0.26619D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.09883D0
	    p(1) = 1.96788D0
	    p(2) = 5.19182D0
	    p(3) = 8.82173D0
	    p(4) = 0.06468D0
	    p(5) = 0.11297D0
	    p(6) = 1.63850D0
	    p(7) = -0.00872D0
	    p(8) = 0.52753D0
	    p(9) = 0.41794D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.14258D0
	    p(1) = 2.00822D0
	    p(2) = 6.23508D0
	    p(3) = 7.81846D0
	    p(4) = 0.07064D0
	    p(5) = -0.05869D0
	    p(6) = 1.24848D0
	    p(7) = -0.01160D0
	    p(8) = 0.48932D0
	    p(9) = 0.22001D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.16184D0
	    p(1) = 2.16963D0
	    p(2) = 7.62378D0
	    p(3) = 7.33369D0
	    p(4) = 0.09197D0
	    p(5) = 0.15692D0
	    p(6) = 1.80734D0
	    p(7) = -0.01561D0
	    p(8) = 0.53224D0
	    p(9) = 0.39357D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.20205D0
	    p(1) = 2.28733D0
	    p(2) = 9.10375D0
	    p(3) = 7.24877D0
	    p(4) = 0.08325D0
	    p(5) = 0.36941D0
	    p(6) = 2.39131D0
	    p(7) = -0.00057D0
	    p(8) = 0.41640D0
	    p(9) = 0.90531D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.95664D0
	    p(1) = 1.11409D0
	    p(2) = 19.00631D0
	    p(3) = 15.97282D0
	    p(4) = 0.15616D0
	    p(5) = 0.40229D0
	    p(6) = 0.85878D0
	    p(7) = -0.03123D0
	    p(8) = 6.75437D0
	    p(9) = -3.83159D0
	  ENDIF

	ENDIF

	off_mKP_fit = -( p(0) * x**p(3) * DEXP(p(1) * x**p(2))
     &                + p(4) * x*DEXP(((x-p(5))/p(6))**2)
     &                + x**0.5D0 * p(7) * DEXP(((x-p(9))/p(8))**2) )

	RETURN
	END
C ***********************************************************************
CCCC   Converts reduced cross sections to structure functions   CCCCC
CCCC   for protons and neutrons                                 CCCCC

      SUBROUTINE SF(w2,q2,F1p,FLp,F2p,F1n,FLn,F2n)
      IMPLICIT none

      real*8 w2,q2,x,sigtp,siglp,sigtn,sigln,f1p,f2p,fLp
      real*8 f1n,f2n,fLn,pi,pi2,alpha,mp,mp2
      Integer i
 

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      x = q2/(q2+w2-mp2)

   
      call rescsp(w2,q2,sigTp,sigLp)
      call rescsn(w2,q2,sigTn,sigLn)

      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      f1n = sigTn/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)

c      fLp = 2*q2/abs(w2-mp2+q2)*sigLp/0.3894e3/pi2/alpha/8.0*
c     &        abs(w2-mp2) 
c      fLn = 2*q2/abs(w2-mp2+q2)*sigLn/0.3894e3/pi2/alpha/8.0*
c     &        abs(w2-mp2) 
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)
      fLn = sigLn*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)


      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)
      f2n = (2.*x*f1n+fLn)/(1.+4.*mp2*x*x/q2)

      return
      end

      
      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis

      data xvalp / 

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
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)
c      call disp(w2,q2,sigTdis,sigLdis)
      if(w2.GT.8.0) then 
c        write(6,*) "resmod: ",w2,q2,sigTp,sigLp
        call disp(w2,q2,sigTdis,sigLdis)
c        write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
        if(w2.LE.10.0) then
         sigTp = (10.0-w2)*sigTp+(w2-8.0)*sigTdis
         sigLp = (10.0-w2)*sigLp+(w2-8.0)*sigLdis
         sigTp = sigTp/2.0
         sigLp = sigLp/2.0
         else
          sigTp = sigTdis
          sigLp = sigLdis
        endif
      endif

      return
      end



      SUBROUTINE RESCSN(w2,q2,sigtn,sigln)
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i

      data xvaln / 
     & 0.12300E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13223E+00,0.21779E+00,0.10504E+00,0.19322E+00,
     & 0.16957E+00,0.38000E+00,0.76000E+01,0.44484E+01,0.22412E+01,
     & 0.18607E+01,0.13787E-02,0.98347E+04,0.18160E+00,0.26770E+01,
     & 0.34272E+00,0.18609E-05,0.47226E+01,0.63725E-06,0.23236E-01,
     & 0.84561E+03,0.29815E+00,0.27082E+01,0.00000E+00,0.10000E+01,
     & 0.10000E+01,0.20000E+01,0.67680E+01,0.56785E+01,0.43946E+00,
     & 0.56784E+01,0.20954E+03,0.18496E+00,0.16769E+01,0.15454E+00,
     & -.83108E+02,0.30574E+01,0.10000E-02,0.93889E+00,-.89302E-02,
     & -.62215E-01,0.19800E+01,0.45000E+00,0.40853E+01,0.14462E+00,
     & 0.10046E+01,0.99612E+00,0.99930E+00,0.99357E+00,0.10192E+01,
     & 0.10000E+01,0.99669E+00,0.10002E+01,0.10000E+01,0.10000E+01,
     & 0.10094E+01,0.10064E+01,0.16270E+03,0.86884E+01,0.57864E+01,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.30131E+02,0.55479E+01,0.93622E+00,0.00000E+00,0.76450E+02,
     & 0.66412E-01,0.85428E+01,0.00000E+00,0.00000E+00,0.20000E+01,
     & 0.10000E+01,0.00000E+00,0.23539E+01,0.10360E-02,0.67452E+00,
     & 0.00000E+00,0.23312E+03,0.35946E+00,0.23000E+01,-.30155E-01,
     & 0.59653E+01,-.10616E+01,0.00000E+00,0.00000E+00,0.13868E+04,
     & 0.29402E+03,0.39655E+00,0.00000E+00,0.26453E+00,0.00000E+00 /


      do i=1,50
        xval1(i) = xvaln(i)
        xvalL(i) = xvaln(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodn(1,w2,q2,xval1,sigtn)
      call resmodn(2,w2,q2,xvalL,sigLn)

      return
      end

      SUBROUTINE RESMODP(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


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
        if(sf.eq.1) sigrsv(i) = sigr(i)
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
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q2+q20)
     &        *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
        enddo

      endif


      sig = sig_res + sig_nr
       
      if(w2.LT.1.159) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 
  

      SUBROUTINE RESMODN(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.939565

c      mp = 0.938272

      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)
      if(w.LT.(mp+mpi)) wdif(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) wdif(2) = 0.0

c      write(6,*) "here"

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
      xpr(2) = 1.+(w2-(mp+2.*mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LT.(mp+mpi)) xpr(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) xpr(2) = 0.0

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
        height(7) = xval(49)*dip

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
        if(sf.eq.1) sigrsv(i) = sigr(i)
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
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
c     &      (1.-xb)**(nr_coef(i,2)+0.0*t )
     &      (1.-xpr(1))**(nr_coef(i,2)+nr_coef(i,3)*log(q2+m0))
c     &        *xpr(1)**(xval(41)+xval(42)*t)
     &        *xb**(xval(41)+xval(42)*log(q2+m0)+
     &                  nr_coef(i,4)*q2*sqrt(q2))
c     &        *(1./(1.+q2/nr_coef(i,3)))**nr_coef(i,4)
        enddo

CCC  next line is test CCC nonres3 directory ******

c        sig_nr = sig_nr + nr_coef(i,3)/(1.+q2/nr_coef(i,4))**2*
c     &               (1.-xpr(1))**2

        

      endif


      sig_res = abs(sig_res)
      sig_nr = abs(sig_nr)

      sig = sig_res + sig_nr
       
      if(w.LT.(mp+mpi)) sig = 0.0

 1000  format(8f12.5)

      RETURN 
      END 

      SUBROUTINE DISP(w2,q2,sigt,sigl)
      IMPLICIT none

      real*8 w2,q2,q2t,x,sigt,sigl,f1,f2,fL,r,dr,r1,r2
      Integer i,modt
      character*1 targ
      real*8 mp2,pi,pi2,alpha,t1,t2
      logical goodfit


      targ = 'P'
      modt = 12
      mp2 = 0.938272
      mp2 = mp2*mp2
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.0359

      x = q2/(w2+q2-mp2)

      call f2allmer(x,q2,f2)
      call f2glober(x,q2,targ,modt,f2)
      call r1998er(x,q2,r,dr,goodfit)
      if(q2.LE.0.15) then
        q2t = 0.15
        call r1998er(x,q2t,r1,dr,goodfit)
        r = r1/0.25*q2
c        write(6,*) r1,r
        call f2allmer(x,q2,f2)
      endif


      f1 = f2/2./x/(r+1.0)*(1.0+4.0*mp2*x*x/q2)
      fL = 2.*x*r*f1    

      sigt = 0.3894e3*f1*pi2*alpha*8.0/abs(w2-mp2)
      sigL = r*sigt

      return
      end
!Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!            Ph.D. Thesis, Stanford University,                                 
!            March 1990.                                                        
!For details see file HELP.DOCUMENT.                                            
                                                                                
!Program contains 145 lines of Fortran code, of 72 characters each, with        
!no subroutines.  Program requires File E.14 as input.                          
                                                                                
                                                                                
      SUBROUTINE F2GLOBER(X,Q2,Target,MODEL,F2)        
                                                                                
! Returns F2 and related quantities from the either the LAMBDA12 model          
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.                   
!                                                                               
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2          
!                                                                               
! Further, program returns uncertainty in F2 based on both statistics           
! (ST) and systematic effects due to our choice of models (SY).                 
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical        
! uncertainty in this slope.                                                    
!                                                                               
! Best model is LAMBDA12.  SY is estimated to be the difference between         
! the two models.                                                               
!                                                                               
! Errors due to overall normalization are not included in program output        
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.              
! Systematic errors due to radiative corrections are shown in Reference         
! to be very kinematic independent, and are everywhere <.5%, and thus,          
! are ignored by this routine (also see documentation to dFRC in file           
! HELP.DOCUMENT).                                                               
!                                                                               
! Coefficients and correlation matrix elements are from File                    
! E.13 F1990.MATRICES, dated 27Jan90.                                           
                                                                                
      IMPLICIT NONE                                                             
      LOGICAL GOODFIT,FIRST                                             
      REAL*8   X, Q2, F2,ST,SY,SLOPE,DSLOPE                               
      REAL*8   XPP,XP,Y,POLY,F2B,POL1,STB,Q,QTH,DQ,F2TH,QUAD,SCB,F2L      
      REAL*8   BINDING,STL                                                
      INTEGER MODEL, I, J, K                                                    
      CHARACTER*1 TARGET                                                        
      REAL*8    B(2,9),L(2,9,9)                      ! OMEGA9 variables         
      REAL*8    C(2,12),M(2,12,12)                   ! LAMBDA12 variable        
      REAL*8 V(9),Z(12),U(12),LIN                                               
                             
      FIRST = .true.                                                   
                                                                                
! Model #9  27 Jan 90.                            !  HYDROGEN Ci            
      DATA (B(1,J),J=1,9)/                        
     >   0.7338659870D0,  11.0245522588D0,   2.6185804129D0,                    
     >   4.0956321483D0,   0.1206495422D0,   1.9714128709D0,                    
     >   3.8893348719D0, -14.0507358314D0,   8.8080576075D0/
!                                                 !  HYDROGEN MijDiDj
      DATA ((L(1,J,K),K=1,9),J=1,9)/              
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,                     
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,                     
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,                     
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,                     
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,                     
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,                     
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,                     
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,                     
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,                     
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,                     
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,                     
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,                     
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,                     
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,                     
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,                     
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,                     
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,                     
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,                     
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,                     
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,                     
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,                     
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,                     
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,                     
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,                     
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,                     
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,                     
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/                     
                                                                                
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci               
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,                     
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,                     
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/                     
                                                                                
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj          
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,                     
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,                     
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,                     
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,                     
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,                     
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,                     
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,                     
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,                     
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,                     
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,                     
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,                     
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,                     
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,                     
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,                     
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,                     
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,                     
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,                     
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,                     
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,                     
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,                     
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,                     
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,                     
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,                     
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,                     
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,                     
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,                     
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/                     
                                                                                
!    MODEL #12:    27Jan90.                                                     
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci                    
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,                     
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,                     
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,                     
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/                     
      DATA ((M(1,J,K),K=1,12),J=1,6)/     !     HYDROGEN MijDiDj               
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,                     
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,                     
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,                     
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,                     
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,                     
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,                     
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,                     
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,                     
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,                     
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,                     
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,                     
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,                     
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,                     
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,                     
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,                     
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,                     
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,                     
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,                     
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,                     
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,                     
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,                     
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,                     
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,                     
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0/                     
        DATA ((M(1,J,K),K=1,12),J=7,12)/             ! hydrogen
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,                     
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,                     
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,                     
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,                     
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,                     
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,                     
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,                     
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,                     
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,                     
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,                     
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,                     
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,                     
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,                     
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,                     
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,                     
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,                     
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,                     
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,                     
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,                     
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,                     
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,                     
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,                     
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,                     
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/                     
                                                                                
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci                   
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,                     
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,                     
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,                     
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/                     
      DATA ((M(2,J,K),K=1,12),J=1,6)/     !     DEUTERIUM MijDiDj              
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,                     
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,                     
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,                     
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,                     
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,                     
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,                     
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,                     
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,                     
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,                     
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,                     
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,                     
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,                     
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,                     
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,                     
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,                     
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,                     
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,                     
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,                     
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,                     
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,                     
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,                     
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,                     
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,                     
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0/                     
        DATA ((M(2,J,K),K=1,12),J=7,12)/             ! deuterium
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,                     
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,                     
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,                     
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,                     
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,                     
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,                     
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,                     
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,                     
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,                     
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,                     
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,                     
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,                     
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,                     
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,                     
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,                     
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,                     
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,                     
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,                     
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,                     
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,                     
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,                     
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,                     
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,                     
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/                       
      !----------------------------------------------------------------         
      !----------------------------------------------------------------         
                                                                                
                                                                                                                          
      i = 1                                                                     
      IF (TARGET.EQ.'D') i = 2                                                  
      BINDING = 1./(1.-EXP(-MIN(20.0,7.7*(1./X+.93828**2/Q2-1.))))               
      IF (i.EQ.1) BINDING = 1.                                                  
                                                                                
      !OMEGA9 MODEL FIRST:                                                      
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                                     
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                                     
           Y    = 1.-XP                                                         
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+              
     >            B(i,9)*Y**7                                                   
           F2B  = X/XPP*BINDING*POLY                                            
          !-----------------------------------------------------------          
          !V(k) is the derivative of F_2 with respect to parameter k.           
           V(1) = -F2B/XPP/(Q2/X+B(i,2))                                        
           V(2) =  F2B/XPP/(Q2/X+B(i,2))**2*(Q2+B(i,1))                         
           POL1 =  3.*B(i,5)*Y**2+4.*B(i,6)*Y**3+5.*B(i,7)*Y**4+                
     >             6.*B(i,8)*Y**5+7.*B(i,9)*Y**6                                
           V(3) = -F2B*POL1/POLY/(Q2/X+B(i,4))                                  
           V(4) =  F2B*POL1/POLY/(Q2/X+B(i,4))**2*(Q2+B(i,3))                   
           DO 10 j = 5,9                                                        
10         V(j) =  F2B/POLY*Y**(j-2)                                            
           STB = 0.                                                             
           DO 11 j = 1,9                                                        
           DO 11 k = 1,9                                                        
11         STB = STB + L(i,j,k)*V(j)*V(k)                                       
           STB = SQRT(STB)*BINDING                                              
                                                                                
      !LAMBDA12 MODEL NEXT:                                                     
           Y    = 1.-X                                                          
           q    = LOG(Q2)                                                       
           qth  = .2+3.2*X                                                      
           dq   = q-qth                                                         
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                          
     >            C(i,4)*Y**6+C(i,5)*Y**7                                       
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                           
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq               
           IF (q.GT.qth) QUAD = 0.                                              
           SCB  = (1.+LIN+QUAD)                                                 
           F2L  = F2th*SCB*BINDING                                              
          !-----------------------------------------------------------          
          !Z(k) is the derivative of F_2 with respect to parameter k.           
           DO 20 j = 1,5                                                        
20         Z(j) = SCB*Y**(j+2)                                                  
           Z(6) = 0.                                                            
           IF (q.LT.qth) Z(6) = F2th*dq**2                                      
           Z(7) = Z(6)*X                                                        
           Z(8) = Z(6)*X**2                                                     
           DO 21 j = 9,12                                                       
21         Z(j) = F2th*X**(j-9)*dq                                              
           STL = 0.                                                             
           DO 22 j = 1,12                                                       
           DO 22 k = 1,12                                                       
22         STL = STL + M(i,j,k)*Z(j)*Z(k)                                       
           STL = SQRT(STL)*BINDING                                              
                                                                                
          !U(k) is the derivative of slope with respect to parameter k.         
           SLOPE= F2th*LIN/dq*BINDING                                           
           DO 30 j = 1,5                                                        
30         U(j) = LIN/dq*Y**(j+2)                                               
           DO 31 j = 6,8                                                        
31         U(j) = 0.                                                            
           DO 32 j = 9,12                                                       
32         U(j) = Z(j)/dq                                                       
           DSLOPE = 0.                                                          
           DO 33 j = 1,12                                                       
           DO 33 k = 1,12                                                       
33         DSLOPE = DSLOPE + M(i,j,k)*U(j)*U(k)                                 
           DSLOPE = SQRT(DSLOPE)                                                
      !----------------------------------------------------------------         
                                                                                
      F2 = 0.                                                                   
      ST = 0.                                                                   
      IF (MODEL.EQ. 9) THEN                                                     
           F2 = F2B                                                             
           ST = STB                                                             
      ELSEIF (MODEL.EQ.12) THEN                                                 
           F2 = F2L                                                             
           ST = STL                                                             
      ELSE                                                                      
           WRITE(*,'('' F1990: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')            
      ENDIF                                                                     
      SY = ABS(F2B-F2L)                                                         
                                                                                
      GOODFIT = .TRUE.                                                          
      !The following cuts define the region of applicability of F1990.          
      !In order they are:                                                       
      !     [radiative corrections convergence criteria in Q2] .and.            
      !     [radiative corrections convergence criteria in x]  .and.            
      !     [stay out of resonance region, W2.ge.3.0]          .and.            
      !     [limitation imposed by maximum beam energy].                        
      IF ((Q2.LT..566).OR.(X.LT..062).OR.(X.LT.Q2/(2.*.93828*21.))              
     >   .OR.(X.GT.1./((3.-.93828**2)/Q2+1.)))     THEN                         
                                                                                
C         WRITE(*,'('' WARNING[F1990]: OUTSIDE RECOMMENDED RANGE.'')')          
          GOODFIT=.FALSE.                                                       
      ENDIF                                                                     
                                                                                
      RETURN                                                                    
      END                                                                       

! Modified by Yongguang Liang to include the R1998 fitting parameters
! 
!File R1990.FORTRN.                                                       
!Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!            Ph.D. Thesis, Stanford University,                                 
!            March 1990.                                                        
!For details see file HELP.DOCUMENT.                                            
                                                                                
!Program contains 135 lines of Fortran code, of 72 characters each, with        
!one subroutine.                                                                
                                                                                
                                                                                
      SUBROUTINE R1998er(X,Q2,R,DR,GOODFIT)                                       
                                                                                
! Model for R, based on a fit to world R measurements. Fit performed by         
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details           
! see Reference.                                                                
!                                                                               
! Three models are used, each model has three free parameters.  The             
! functional forms of the models are phenomenological and somewhat              
! contrived.  Each model fits the data very well, and the average of            
! the fits is returned.  The standard deviation of the fit values is            
! used to estimate the systematic uncertainty due to model dependence.          
!                                                                               
! Statistical uncertainties due to fluctuations in measured values have         
! have been studied extensively.  A parametrization of the statistical          
! uncertainty of R1990 is presented in FUNCTION DR1990.                         
!                                                                               
! The three model forms are given by:                                           
!                                                                               
!     R_A = A(1)/LOG(Q2/.04)*FAC + A(2)/[Q2ã4+A(3)ã4]ã.25 ;                     
!     R_B = B(1)/LOG(Q2/.04)*FAC + B(2)/Q2 + B(3)/(Q2**2+.3**2) ;               
!     R_C = C(1)/LOG(Q2/.04)*FAC + C(2)/[(Q2-Q2thr)ã2+C(3)ã2]ã.5 ,              
!                               ...where Q2thr = 5(1-X)ã5 ;                     
!           where FAC = 1+12[Q2/(1+Q2)][.125ã2/(.125ã2+xã2)] gives the          
!           x-dependence of the logarithmic part in order to match Rqcd         
!           at high Q2.                                                         
!                                                                               
! Each model fits very well.  As each model has its own strong points           
! and drawbacks, R1990 returns the average of the models.  The                  
! chisquare for each fit (124 degrees of freedom) are:                          
!     R_A: 110,    R_B: 110,    R_C: 114,    R1990(=avg): 108                   
!                                                                               
! This subroutine returns reasonable values for R for all x and for all         
! Q2 greater than or equal to .3 GeV.                                           
!                                                                               
! The uncertainty in R originates in three sources:                             
!                                                                               
!     D1 = uncertainty in R due to statistical fluctuations of the data         
!          and is parameterized in FUNCTION DR1990, for details see             
!          Reference.                                                           
!                                                                               
!     D2 = uncertainty in R due to possible model dependence, approxi-          
!          mated by the variance between the models.                            
!                                                                               
!     D3 = uncertainty in R due to possible epsilon dependent errors            
!          in the radiative corrections, taken to be +/- .025.  See             
!          theses (mine or Dasu's) for details.                                 
!                                                                               
! and the total error is returned by the program:                               
!                                                                               
!     DR = is the total uncertainty in R, DR = sqrt(D1ã2+D2ã2+D3ã2).            
!          DR is my best estimate of how well we have measured R.  At           
!          high Q2, where R is small, DR is typically larger than R.  If        
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,          
!          you might wish to assume DR = 0 at very high Q2.                     
!                                                                               
! NOTE:    In many applications, for example the extraction of F2 from          
!          measured cross section, you do not want the full error in R          
!          given by DR.  Rather, you will want to use only the D1 and D2        
!          contributions, and the D3 contribution from radiative                
!          corrections propogates complexely into F2.  For more informa-        
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or             
!          for explicite detail, see Reference.                                 
!                                                                               
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      IMPLICIT NONE                                                             
      REAL*8 FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1990,X,Q2,W2                
      REAL*8 A(3), B(3), C(3), MP,MP2
      LOGICAL GOODFIT                                                           
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      DATA A/ .06723, .46714, 1.89794 /                                  
      DATA B/ .06347, .57468, -.35342 /                                   
      DATA C/ .05992, .50885, 2.10807 /

      MP = .9382727
      MP2 = MP*MP                                                                        

      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))                       
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!           
c      Q2thr = 5.*(1.-X)**5                                                      
                                                                                
c      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))                        
c      R_B   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)                          
c      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)                      


       W2 = Q2/X-Q2+MP2

c       IF (W2 .GT. 3.9 ) THEN 
       q2thr = 12.3708*x-43.1043*x**2  !! r1998
     & +41.7415*x**3

      r_a = 0.0485*rlog
     &+0.5470/sqrt(sqrt(q2**4+2.0621**4))
     &*(1.-0.3804*x+0.5090*x**2)*x**(-0.0285) 
      r_b = 0.0481*rlog
     &+(0.6114/q2-0.3509/(q2**2+.3**2))
     &*(1.-0.4611*x+0.7172*x**2)*x**(-0.0317)                    
      r_c = 0.0577*rlog
     & + 0.4644/sqrt((q2-q2thr)**2+1.8288**2)
c      ELSE
c       q2thr = -5.8119*x+0.9259*x**2  !! r1990
c     & +9.4095*x**3

c      r_a = 0.2839*rlog
c     &+0.0789/sqrt(sqrt(q2**4+0.0**4))
c     &*(1.+32.2507*x-28.5938*x**2)*x**(2.8009) 
c      r_b = 0.2518*rlog
c     &+(0.0110/q2-0.0069/(q2**2+.3**2))
c     &*(1.+57.0509*x-1.4183*x**2)*x**(-0.0934)                    
c      r_c = 0.0621*rlog
c     & + 1.4431/sqrt((q2-q2thr)**2+7.1946**2)
c       write(*,*) r_a, r_b, r_c 
c      ENDIF
      
      R     = (R_A+R_B+R_C)/3.                                                  
                                                                                
      D1    = DR1990(X,Q2)                                                      
      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)                       
      D3    = .023*(1.+.5*R)                                                    
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                               
      DR    = SQRT(D1**2+D2**2+D3**2)                                           
                                                                                
      GOODFIT = .TRUE.                                                          
      IF (Q2.LT..3) GOODFIT = .FALSE.                                           
      RETURN                                                                    
      END                                                                       
                                                                                
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
                                                                                
      FUNCTION DR1990(X,Q2)                                                     
                                                                                
! Parameterizes the uncertainty in R1990 due to the statistical                 
! fluctuations in the data.  Values reflect an average of the R-values          
! about a neighborhood of the specific (x,Q2) value.  That neighborhood         
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see               
! Reference.                                                                    
!                                                                               
! This subroutine is accurate over all (x,Q2), not only the SLAC deep           
! inelastic range.  Where there is no data, for example in the resonance        
! region, it returns a realistic uncertainty, extrapolated from the deep        
! inelastic region (suitably enlarged).  We similarly estimate the              
! uncertainty at very large Q2 by extrapolating from the highest Q2             
! measurments.  For extremely large Q2, R is expected to fall to zero,          
! so the uncertainty in R should not continue to grow.  For this reason         
! DR1990 uses the value at 64 GeV for all larger Q2.                            
!                                                                               
! XHIGH accounts for the rapidly diminishing statistical accuracy for           
! x>.8, and does not contribute for smaller x.                                  
                                                                                
                                                                                
      IMPLICIT NONE                                                             
      REAL*8 U(10,10),DR1990,QMAX,Q,S,A,XLOW,XHIGH,X,Q2                           
                                                                                
                                                                                
      QMAX = 64.                                                                
                                                                                
      Q = MIN(Q2,QMAX)                                                          
      S = .006+.03*X**2                                                         
      A = MAX(.05,8.33*X-.66)                                                   
                                                                                
      XLOW  = .020+ABS(S*LOG(Q/A))                                              
      XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)                          
                                                                                
      DR1990 = SQRT(XLOW**2+XHIGH**2)                                           
      RETURN                                                                    
      END                                                                       
!                                                                               
!End of file R1990.FORTRN.  135 Fortran lines.                                  





c       allm97, NMC published measured points Q2>0.75 GeV2
c       for values Q<1 use data of E665!
c       parameterization of F2 , according to
c       H.Abramowicz and A.Levy, hep-ph/9712415
c
c       3*10-6 < x  < 0.85, W2>3GeV2
c       0.   < Q2 < 5000 GeV2, dof=0.97
c
 
      SUBROUTINE f2allmer(x,q2,f2a)

      IMPLICIT NONE
 
      REAL*8 x,q2,M22,f2a
      REAL*8 SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      REAL*8 S11,A11,B11,M12,S21,A21,B21,M02,LAM2,Q02
      REAL*8 S12,S13,A12,A13,B12,B13,S22,S23,A22,A23
      REAL*8 B22,B23,w2,w,z
      REAL*8 ALFA,XMP2
C  POMERON
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )         
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      
C
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
        
C
      IF(Q2.EQ.0.) THEN
       S=0.
       Z=1.

C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11
       BP=B11
       SP=S11
       F2P=SP*XP**AP
C                                               
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21
       BR=B21
       SR=S21
       F2R=SR*XR**AR
C
      ELSE
       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))
       Z=1.-X   
C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)
       BP=B11+B12*S**B13
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)
       F2P=SP*XP**AP*Z**BP
C
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21+A22*S**A23
       BR=B21+B22*S**B23
       SR=S21+S22*S**S23
       F2R=SR*XR**AR*Z**BR
 
C
      ENDIF                                     
      
 
      f2a = q2/(q2+m02)*(F2P+F2R)
 
 
      RETURN
      END                                  
