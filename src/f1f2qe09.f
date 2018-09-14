      SUBROUTINE F1F2QE09(Z, A, qsq, wsq, F1, F2)
c
C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS
c for A>2 uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c for A=2 uses pre-integrated Paris wave function (see ~bosted/smear.f)
c coded by P. Bosted August to October, 2006
c
c input: Z, A  (real*8) Z and A of nucleus (shoud be 2.0D0 for deueron)
c        Qsq (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        Wsq (real*8) is invarinat mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus
c
c Note: Deuteron agrees well with Laget (see ~bosted/eg1b/laget.f) for
c a) Q2<1 gev**2 and dsig > 1% of peak: doesnt describe tail at high W
c b) Q2>1 gev**2 on wings of q.e. peak. But, this model is up
c    to 50% too big at top of q.e. peak. BUT, F2 DOES agree very
c    nicely with Osipenko et al data from CLAS, up to 5 GeV**2

      IMPLICIT NONE     
      REAL*8 P(0:23)
      COMMON/PARCORR/P
      REAL*8 Z, A, avgN, F1, F2, wsq, qsq
      REAL*8 amp/0.93828/, amd/1.8756/
      REAL*8 PAULI_SUP1, PAULI_SUP2
      REAL*8 GEP, GEN, GMP, GMN, Q, Q3, Q4
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      real*8 pz, nu, dpz, pznom, pzmin
      real*8 qv, TAU, W1, W2, FY, dwmin, w2p
      real kappa, lam, lamp, taup, squigglef, psi, psip, nuL, nut
      real kf, es, GM2bar, GE2bar, W1bar, W2bar, Delta, GL, GT
      real*8 mcor,ecor
      integer IA, izz, izzmin, izp, izznom, izdif
      
      real*8 a1,a2,a3,b1,b2,b3,b4,b5,gd

      real*8 xvalc(40) /
     & 0.47352E+01,0.11441E+01,0.24509E+00,0.53496E+01,0.86248E+00,
     & 0.13058E+00,0.13429E+01,0.17119E+01,0.19525E+00,-.16900E+01,
     & 0.98573E+00,0.97944E+00,0.10449E+01,0.99495E+00,0.98709E+00,
     & 0.10000E+01,0.97853E+00,0.10066E+01,0.98571E+00,0.99956E+00,
     & 0.98897E+00,0.99984E+00,0.10006E+01,0.10000E+01,0.10000E+01,
     & 0.10800E+01,0.77077E+00,0.57666E+00,0.28701E+00,0.64186E+02,
     & 0.83304E+01,0.16978E-01,0.10841E+00,0.15743E+00,-.11446E+00,
     & 0.17658E+01,0.10000E+01,0.52408E+02,0.14887E+00,0.86128E+01 /

c Look up tables for deuteron case
       real*8 fyd(200)/
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
       real*8 avp2(200)/
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

c     Peter Bosted's correction params

       real*8 pb(20)/ 0.1023E+02, 0.1052E+01, 0.2485E-01, 0.1455E+01,
     >      0.5650E+01,-0.2889E+00, 0.4943E-01,-0.8183E-01,
     >     -0.7495E+00, 0.8426E+00,-0.2829E+01, 0.1607E+01,
     >      0.1733E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,
     >      0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00/ 
       real*8 y,R

! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN

! some kinematic factors. Return if nu or qsq is negative
      Nu = (wsq - amp**2 + qsq) / 2. / amp
      if(nu .le. 0.0 .or. qsq .lt. 0.) return
      TAU   = QSQ / 4.0 / amp**2                                        
      qv = sqrt(nu**2 + qsq)

! Bosted fit for nucleon form factors Phys. Rev. C 51, p. 409 (1995)
      Q = sqrt(QSQ)
      Q3 = QSQ * Q
      Q4 = QSQ**2
      GEP = 1./  (1. + 0.14 * Q + 3.01 * QSQ + 0.02 * Q3 + 
     >  1.20 * Q4 + 0.32 * Q**5)
      GMP = RMUP * GEP
      GMN = RMUN / (1.- 1.74 * Q + 9.29 * QSQ - 7.63 * Q3 + 
     >  4.63 * Q4)
      GEN = 1.25 * RMUN * TAU / (1. + 18.3 * TAU) / 
     >  (1. + QSQ / 0.71)**2

CCCCCCCCCCCCCCCCCCC  New FormFactors  CCCCCCCCCCCCCCCCCCCCCCCCC

c        GMp =2.792782/(1.-0.29940*q+5.65622*qsq-5.57350*q*qsq+      !  with Hall A data  !
c     &         6.31574*qsq**2-1.77642*q*qsq**2+0.30100*qsq**3)

c        GEp = 1./(1.-0.31453*q+7.15389*qsq-14.42166*q*qsq+
c     &         23.33935*qsq**2-14.67906*q*qsq**2+4.03066*qsq**3)

c        GMp = 2.792782     !!!    
c     &       * (1.D0 -1.465D0*tau + 1.260D0*tau**2 + 0.262D0*tau**3)
c     &       / (1.D0 + 9.627D0*tau + 0.0D0*tau**2 + 0.0D0*tau**3
c     &               + 11.179D0*tau**4 + 13.245D0*tau**5)


        a1 = -1.651D0
        a2 = 1.287D0
        a3 = -0.185D0
        b1 = 9.531D0
        b2 = 0.591D0
        b3 = 0.D0
        b4 = 0.D0
        b5 = 4.994D0

        GEp = (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)


        a1 = -2.151D0
        a2 = 4.261D0
        a3 = 0.159D0
        b1 = 8.647D0
        b2 = 0.001D0
        b3 = 5.245D0
        b4 = 82.817D0
        b5 = 14.191D0

        GMp = RMUP
     &      * (1.D0 + a1*tau + a2*tau**2 + a3*tau**3)
     &      / (1.D0 + b1*tau + b2*tau**2 + b3*tau**3
     &              + b4*tau**4 + b5*tau**5)

        GD = (1./(1 + qsq/0.71))**2
        GMn = -1.913148           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &  /(1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)

*        GMn = -1.913148*GD

        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        mcor = 1.0-0.09*qsq/(0.1+qsq)
        ecor = 1.0+0.15*qsq/(0.1+qsq)
        GMn = mcor*GMn
        GEn = ecor*GEn


! Get kf and Es from superscaling from Sick, Donnelly, Maieron,
c nucl-th/0109032
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
c         if(iA.gt.7) Es=0.015
c      if(IA.gt.7) kf = 0.228
c      if(iA.gt.7) Es= 0.015

c Test MEC  10/13
         if(iA.gt.7) Es=0.015
         if(IA.gt.7) kf=0.228
c

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


! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

! structure functions with off shell factors
      kappa = qv / 2. / amp
      lam = nu / 2. / amp
      lamp = lam - Es / 2. / amp
      taup = kappa**2 - lamp**2
      squigglef = sqrt(1. + (kf/amp)**2) -1.
! Very close to treshold, could have a problem
      if(1.+lamp.le.0.) return
      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - tau ) / sqrt(squigglef) /
     >  sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))

      psip = (lamp - taup) / sqrt(squigglef) / 
     >  sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))

      nuL = (tau / kappa**2)**2

c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 

      GM2bar = Pauli_sup1 * (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = Pauli_sup2 * (Z * GEP**2 + avgN * GEN**2) 
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

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
c      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
c     >   (1. + exp(-2.4291 * psip)) / kf

CCMEC - testing below  CCCC

      FY = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >   (1. + exp(xvalc(10) * psip)) / kf


! Use PWIA and Paris W.F. for deuteron to get better FY
      if(IA.eq.2) then
! value assuming average p2=0.
        pz = (qsq - 2. * amp * nu ) / 2. / qv
        izz = int((pz + 1.0) / 0.01) + 1
        izz = min(200,max(1,izz))
        izznom = izz
! ignoring energy term, estimate change in pz to compensate
! for avp2 term
        dpz = avp2(izznom) / 2. / qv
        izdif = dpz * 150. 
        dwmin=1.E6
        izzmin=0
        do izp = izznom, min(200, max(1, izznom + izdif))
          pz = -1. + 0.01 * (izp-0.5)
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izp)))**2 - 
c    >      qv**2 + 2. * qv * pz - avp2(izp)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izp)
c if passed first minimum, quit looking so don't find second one
          if(abs(w2p - amp**2).gt.dwmin) goto 11
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            izzmin = izp
          endif
        enddo
 11     izz = min(199,max(2,izzmin))
! search for minimum in 1/10th bins locally
        pznom = -1. + 0.01 * (izz-0.5)
        dwmin=1.E6
        do izp = 1,19
          pz = pznom - 0.01 + 0.001 * izp
c *** this version gives worse agreement with laget than
c         w2p = (amd + nu - sqrt(amp**2 + avp2(izz)))**2 - 
c   >      qv**2 + 2. * qv * pz - avp2(izz)
c this version!
          w2p = (amd + nu - amp )**2 - 
     >      qv**2 + 2. * qv * pz - avp2(izz)
          if(abs(w2p - amp**2).lt.dwmin) then
            dwmin = abs(w2p - amp**2)
            pzmin = pz
          endif
        enddo
        if(dwmin.ge.1.e6.or.abs(pznom-pzmin).gt.0.01) 
     >     write(6,'(1x,''error in dwmin,pzmin'',3i4,6f7.3)')
     >     izznom,izzmin,izz,qsq,wsq,w2p,dwmin/1.e6,pzmin,pznom
        if(pzmin.lt.pznom) then
          fy = fyd(izz) - (fyd(izz-1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        else
          fy = fyd(izz) + (fyd(izz+1) - fyd(izz)) * 
     >      (pzmin - pznom) / 0.01
        endif
      endif

c final results
      F2 = nu * FY * (nuL * GL + nuT * GT)
      F1 = amp * FY * GT / 2.

      if(F1.LT.0.0) F1 = 0.
      if(nu.gt.0. .and.f1.gt.0.) then
        R = (F2 / nu) / (F1 / amp) * (1. + nu**2 / qsq) - 1.0
      else
        r = 0.4/qsq
      endif


c apply correction factors
      if(A.gt.2) then
         y = (wsq -amp**2) / qv
c         F1 = F1 * (1. + pb(8) + pb(9) * y +
c     >        pb(10)*y**2 + pb(11)*y**3 + pb(12)*y**4 )
c         R = R * (1. + pb(13))
c         F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)

cc correction to correction Vahe
         if(wsq.gt.0.0) then

c            F1=F1*(1.0+P(7)+P(8)*y+P(9)*y**2 +P(10)*y**3 +P(11)*y**4)
c            R = R * ( 1.0 + P(12) )
            F2 = nu * F1/amp * (1. + R) / (1. + nu**2/qsq)
            if(F1.LT.0.0) F1=0.0

         endif
      endif

      return
      end

      BLOCK DATA CORRECTION
      REAL*8 P(0:23)
      COMMON/PARCORR/P
c     DATA P/
c    c       5.1141e-02,   9.9343e-01,   5.3317e-02,   1.3949e+00, 
c    c       5.9993e+00,  -2.9393e-01,   9.9316e-02,   1.0935e-02, 
c    c       3.7697e-01,   1.3542e+00,  -7.4618e+00,   4.3540e+00, 
c    c      -3.7911e-01,   3.5105e-01,  -1.8903e+00,   4.9139e+00, 
c    c      -5.9923e+00,   2.5021e+00,   1.9943e+01,  -5.5879e-02, 
c    c       3.1914e-01,  -1.8657e-01,   3.2746e+01,   4.9801e-03/
       DATA P/
     c       5.1377e-03,   9.8071e-01,   4.6379e-02,   1.6433e+00,
     c       6.9826e+00,  -2.2655e-01,   1.1095e-01,   2.7945e-02,
     c       4.0643e-01,   1.6076e+00,  -7.5460e+00,   4.4418e+00,
     c      -3.7464e-01,   1.0414e-01,  -2.6852e-01,   9.6653e-01,
     c      -1.9055e+00,   9.8965e-01,   2.0613e+02,  -4.5536e-02,
     c       2.4902e-01,  -1.3728e-01,   2.9201e+01,   4.9280e-03/
       END
