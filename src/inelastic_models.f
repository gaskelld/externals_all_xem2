C DJG: This file contains various older, inelastic models (F2, cross section, R) found in "externals_all"
C as of 9/13/2018.
C Most of these routines are kept for historical purposes. Not used so much these days.
C 

!-----------------------------------------------------------------------
      SUBROUTINE F2_JOIN(QQ,X,amuM,iZ,INEL_MODEL,W1,W2,iA) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon by joining the resonance
!  model with the DIS model
! amuM=atomic weight
! INEL_MODEL =rrdd where rr is resonance model number and dd is DIS model 
!             It must be >=100
! DIS_MODEL= 1 ineft
!            2 f2nmc
!            3 f2nmc95
!            4 SMC
!            5 E665
!            6 f2allm (HERA fit, obtained from Antje Burrel
!            9 f2glob model 9
!           12 f2glob model 12
!           33 rock_res9
! RES_MODEL  1 ineft
!            2 H2Model
!            3 rock_res9
!            5 E665
!            
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
! EMC effect corrected for neutron excess 8/19/98
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,ERR_LO,ERR_HI,NU,WSQ,f2allm
!      real fnp_nmc,corr
      real w1model,w2model,w1in,w2in,fnp_nmc
      INTEGER IZ,iA
      REAL MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC_N,FRAC,amuM,W1_DIS,W2_DIS,DSLOPE,ST,SLOPE,SY
      real f2nf2p,emcfacp,ROCK_RES9
      REAL*8 W18,W28,QQ8,W8,amuM8,F28,RDUM
      INTEGER INEL_MODEL,IHD,RES_MODEL,DIS_MODEL
      LOGICAL GD
      CHARACTER*1 HD(2)/'H','D'/
!** changed from 0.3 to 0.0 to avoid discontinuities. H2model
!** is well enough behaved that this should be fine
      REAL Q_RES_LIMIT/0.0/  
      logical ReallyJoin
      common/testing/prttst,usegrd
      logical prttst,usegrd

      RES_MODEL =INEL_MODEL/100
      DIS_MODEL =INEL_MODEL -RES_MODEL*100

      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      QQ8=QQ
      W8= SQRT(max(0.,WSQ))
      amuM8 =amuM 
      IHD =MIN(IA,2)   !no more than deuteron

      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC_N(X,REAL(IA),REAL(IZ),GD)) !with neutron excess 8/19/98
!*** If using F2ALL, use R1990 because this used to get these F2 values
!*** changed 1/19/02
      IF(DIS_MODEL.EQ.6) THEN
        CALL R1990F(X,QQ,R,DR)
      ELSE
        CALL R1998(X,QQ,R,DR,GD)
      ENDIF

!* If A>2 and doing Anje fit, only use F2ALL
      ReallyJoin=.true.
      if(IA.gt.2.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 
!* If RES_MODEL and DIS_Model both 6, don't join
      if(RES_MODEL.eq.6.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 

!** Don't join if using rock_res9: should be good everywhere
      if(DIS_MODEL.EQ.33) ReallyJoin=.false.

      IF(WSQ.GT.3.5.or.(.not.ReallyJoin)) THEN        !DIS
	IF(DIS_MODEL.EQ.1) THEN   !Bodek fit  REAL*8      
          CALL INEFT(QQ8,W8,W18,W28) !already has EMC correction
          W1=W18
          W2=W28
        ENDIF
        IF((DIS_MODEL.GE.2.AND.DIS_MODEL.LE.6).or.
     >     DIS_MODEL.EQ.33 ) THEN
         IF(DIS_MODEL.EQ.2) CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
         IF(DIS_MODEL.EQ.3) CALL F2NMC_NEW(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.4) CALL F2SMC98(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.5)  THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          F2 = F28
         ENDIF
         IF(DIS_MODEL.EQ.33) F2=ROCK_RES9(ihd,x,qq)
         W2 = F2/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF
!** Changed 1/19/02 to use Anje's n/p fit and nuclear dep.
!** Changed 7/20/06 back to NMC for n/p (works MUCH better!)
        IF(DIS_MODEL.EQ.6) THEN  !(HERA fit, obtained from Antje Burrel
          F2 = F2ALLM(X,QQ)   ! proton fit
          IF(iHD.EQ.2)  F2 = (FNP_NMC(X,QQ)+1.)/2. *F2
! This is what Antje uses:
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))
c          IF(IHD.EQ.2)  F2 = (F2NF2P+1.)/2. *F2
!* For C, Al, Cu, Au, override EMCFAC
          emcfacp=emcfac
          if(ia.eq.12) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** use carbon for Al since Antje's seems to low at high x
          if(ia.eq.27) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** turn off Antje's Al corr., 
!         if(ia.eq.27) emcfacp=(0.825-0.46*x-0.19*exp(-21.8*x)+             !!!ALUMINUM  !!!
!    >                    (0.34*x**(-4.91))*x**5.0)
          if(ia.eq.65) emcfacp=(1.026-0.56*x-0.34*exp(-45.7*x)+             !!!COPPER   !!!
     >                    (0.26*x**(-4.41))*x**5.0) 
          if(ia.eq.197) emcfacp=(0.970-1.433*x-0.334*exp(-54.53*x)+          !!!GOLD    !!!
     >                  1.074*x**0.711)
!**** changed to use Rock emc, not Antje!
!***          W2 = F2/NU *EMCFACP
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
          if(prttst) write(*,'(1x,''W2,F2,R,emcfacp='',4f10.4)')
     >      W2,F2,R,emcfacp  
        ENDIF

        IF((DIS_MODEL.EQ.9).OR.(DIS_MODEL.EQ.12)) THEN
          CALL F2GLOB(X,QQ,HD(IHD),DIS_MODEL,F2,ST,SY,SLOPE,DSLOPE,GD)
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R)
        ENDIF
        if(prttst) write(*,'(1x,''model,x,q2,F2,r='',i3,5f10.4)')
     >    DIS_MODEL,x,qq,f2,r,emcfac  
      ENDIF

      IF(WSQ.LT.4.3.and.ReallyJoin) THEN !Resonance
        W1_DIS =W1
        W2_DIS =W2
        IF(RES_MODEL.EQ.1) THEN  !old Bodek resonance fit: EMC corrected
          CALL INEFT(QQ8,W8,W18,W28)
          W1=W18
          W2=W28
        ELSEIF(RES_MODEL.EQ.5) THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          W2 = F28/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF 
        IF(RES_MODEL.EQ.2) THEN !Stuart Fit to Resonance, valid forQ2>.75
          CALL INEFT(QQ8,W8,W18,W28) !already had EMC, neutron corrected
          W1=W18
          W2=W28
          w1in = w1
          w2in = w2
c old way          IF((QQ.LT.10.).AND.(QQ.GE.Q_RES_LIMIT).AND.iA.lt.3) THEN
          IF(QQ.GE.0.3.AND.iA.lt.3) THEN
            if(ia.eq.1) CALL H2MODEL(QQ,WSQ,W1model,W2model)
            if(ia.eq.2) CALL D2MODEL_IOANA(QQ,WSQ,W1model,W2model)
            if(prttst) WRITE(*,'('' H2/D2Model,q2,wsq,w1,w2='',
     >        5f8.3)') QQ,WSQ,W1model,W2model,EMCFAC 
            if(qq.gt.0.5.and.qq.lt.15.) then
              w1=w1model
              w2=w2model
            endif
            if(qq.ge.0.3.and.qq.le.0.5) then
              frac=(qq-0.3)/(0.5-0.3)
              w1 = w1model * frac + w1in *(1.0-frac)
              w2 = w2model * frac + w2in *(1.0-frac)
            endif
            if(qq.ge.10.0.and.qq.le.15.) then
              frac=(qq - 10.0)/(15.0 - 10.0)
              w1 = w1in * frac + w1model *(1.0-frac)
              w2 = w2in * frac + w2model *(1.0-frac)
            endif
          ENDIF
        ENDIF 

        IF(WSQ.GT.3.5) THEN   !interpolate
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           if(prttst) WRITE(*,'(1x,''joining w2'',3f10.3)') 
     >       w2_dis,w2,frac
           W1 = W1_DIS*FRAC + W1*(1.0 - FRAC)
           W2 = W2_DIS*FRAC + W2*(1.0 - FRAC)
        ENDIF
      ENDIF
!**** If using Antje model, then also use Blok correction factors
!* taken out 1/19/02
!      if(DIS_MODEL.EQ.6) THEN
!        if(ihd.eq.1) call f2corr1h(x,corr)
!        if(ihd.eq.2) call f2corr2h(x,corr)
!        w1=w1*corr
!        w2=w2*corr
!      endif
      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE NMCDIS(QQ,X,W1,W2,amuM) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon. 
!  NMC fit for DIS and Old Bodek fit for Resonance.
! IA=atomic weight
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,NU,WSQ, MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC,FRAC,amuM
      REAL*8 W18,W28,QQ8,W8,amuM8
      INTEGER IA,IHD
      LOGICAL GOOD


      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      IF(WSQ.LT.1.16) RETURN
      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC(X,amuM,GOOD)) !isoscaler nucleus
      CALL R1998(X,QQ,R,DR,GOOD)
      IA = amuM+.5  !make into integer
      amuM8 =amuM
      IF(WSQ.GT.3.5) THEN
        
        IHD =MIN(IA,2)   !no more than deuteron
        CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
        W2 = F2/NU *EMCFAC
        W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
      ENDIF

      IF(WSQ.LT.4.3) THEN
         QQ8=QQ
         W8= SQRT(max(0.,WSQ))
         IF(WSQ.LE.3.5) THEN
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit: EMC corrected
           W1=W18
           W2=W28
         ELSE  !interpolate
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit       
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + W18*(1.0 - FRAC)
           W2 = W2*FRAC + W28*(1.0 - FRAC)
         ENDIF
      ENDIF
      F2 = NU*W2
      RETURN
      END

!---------------------------------------------------------------------
       SUBROUTINE F2DGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(deuteron) for E665                          C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2DHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2DHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2DHIW(QQ,NU)

c.....low W parametrization for deuteron
      LOWWF2 = F2LOWW(2,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END

      SUBROUTINE F2PGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(proton) for E665                            C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2PHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2PHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2PHIW(QQ,NU)

c.....low W parametrization for proton
      LOWWF2 = F2LOWW(1,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END
      DOUBLE PRECISION FUNCTION F2HI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** hydrogen BCDMS-like function
      DATA DPEMC /-0.1011D0,2.562D0,0.4121D0,-0.518D0,5.957D0,
     &           -10.197D0,4.695D0,
     &             0.364D0,-2.764D0,0.015D0,0.0186D0,-1.179D0,
     &             8.24D0,-36.36D0,47.76D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2HI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION F2DI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** deuterium BCDMS-like function
      DATA DPEMC /-0.0996D0,2.489D0,0.4684D0,-1.924D0,8.159D0,
     &            -10.893D0,4.535D0,
     &              0.252D0,-2.713D0,0.0254D0,0.0299D0,-1.221D0,
     &              7.5D0,-30.49D0,40.23D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2DI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)

C
      RETURN
      END

      REAL FUNCTION F2DHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(deuteron) for W>5                                *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2DDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2DI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2DI15,F2DDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2DI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2DDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2DHIW = F2

      RETURN
      END
      REAL FUNCTION F2PHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) for W>5                                  *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2PDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2HI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2HI15,F2PDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2HI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2PDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2PHIW = F2

      RETURN
      END
      REAL FUNCTION F2DDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(deuteron)                                        *
* Uses n/p parametrization of E665 and model of Donnachie and*
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
* currently sets n/p to 1                                    *
*                                                            *
* 21may94 AVK                                                *
* set n/p = 0.94                                             *
* private communication of E665 measurement at low x and Q2  *
* from P. G. Spentzouris                                     *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)
      F2 = F2*(1.0+0.94)/2.0

      IF (F2.LT.0.0) F2 = 0.0
      F2DDL = F2

      RETURN
      END
      REAL FUNCTION F2PDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)

      IF (F2.LT.0.0) F2 = 0.0
      F2PDL = F2

      RETURN
      END
      REAL FUNCTION F2DOLA(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
*                                                            *
* The form of the function is a sum of two Regge powers, and *
* is fit to photoproduction data over a wide range of energy *
* and to NMC data for Q2<10 Gev2                             *
* Only the simple function represented by Eqns 4a and 4b are *
* used here, the refinements discussed in the paper are      *
* omitted. F2(neutron) is also omitted.                      *
*                                                            *
* units of Q2ext and Nuext are GeV^2, GeV                    *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,A,SMALLA,B,SMALLB,POWER1,POWER2,TWOM,F2
      PARAMETER (A=0.324,B=0.098,SMALLA=0.562,SMALLB=0.01113)
      PARAMETER (POWER1=0.0808,POWER2=-0.4525)
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      REAL TERM1,TERM2,TWOMNU
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

      TWOMNU = TWOM*NU

      TERM1= A * TWOMNU**POWER1 / (Q2+SMALLA)**(POWER1+1.0)

      TERM2= B * TWOMNU**POWER2 / (Q2+SMALLB)**(POWER2+1.0)

      F2 = Q2 * (TERM1 + TERM2)

      IF (F2.LT.0.0) F2 = 0.0
      F2DOLA = F2

      RETURN
      END
      REAL FUNCTION F2LOWW(ITARGE,Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*   A V Kotwal 30dec93                                       *
*                                                            *
*Return F2 for 1.1<W<5 over a large range of Q2 down to Q2=0 *
*                                                            *
* Based on parametrizations of SLAC,DESY and Daresbury data. *
* Calls two other routines for the resonance region 1.1<W<2  *
* and the inelastic region 2<W<5                             *
* Assume proton and deuteron same for resonance region       *
*                                                            *
* ITARG = 1 for H2, 2 for D2                                 *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
*                                                            *
*     PGS Feb 20, 94                                         *
* Proton and neutron are not the same, see Close chapter 7.2 *
* Modify deuteron resonances accordingly.                    *
* Simple scaling of the D resonance parametrization          *
* according to data (flat line fit in 8 points) by 0.88      *
* Data from J. Franz et al Z.Phys.C 96 P. and F. 10,         *
* 105-116(1981)                                              *
* fit a flat line in the H/D plots p. 110                    *
**************************************************************
      REAL M2,TWOM,F2,Q2,NU,W2,W,WMIN,WINTER,WMAX,NUMAX,WFUNC
      PARAMETER (M2=0.8803,TWOM=1.87654) ! protonmass**2, 2*protonmass
      PARAMETER (WMIN=1.12,WINTER=1.95,WMAX=5.0)
      INTEGER ITARGE,ITARG
      REAL F2LONU,F2PRES,LONUF2,PRESF2,WDIFF
      EXTERNAL F2LONU,F2PRES
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      ITARG = ITARGE

      W2 = M2 + TWOM*NU - Q2
      IF (W2.LT.M2)  W2=M2
      W  = SQRT( W2 )

      IF (NU.LT.0.001) NU = 0.001
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (ITARG.GT.2) ITARG = 2
      IF (ITARG.LT.1) ITARG = 1

      IF (W.LT.WMIN) THEN
c.......force to zero smoothly as W->M
        PRESF2 = F2PRES(Q2,WMIN)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        F2 = PRESF2 * ((W2-M2)/(WMIN**2-M2))**2
      ELSEIF (W.LE.WMAX) THEN
c.....inelastic region, joined smoothly to resonance region
        WDIFF = (W-WINTER)/0.03
        IF (WDIFF.GT.10.0) THEN
          WFUNC = 0.0
        ELSEIF (WDIFF.LT.-10.0) THEN
          WFUNC = 1.0
        ELSE
          WFUNC = 1.0/(1.0+EXP(WDIFF))
        ENDIF
        PRESF2 = F2PRES(Q2,W)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        LONUF2 = F2LONU(ITARG,Q2,NU)
        F2=WFUNC*PRESF2+(1.0-WFUNC)*LONUF2
      ELSE
C.......force to zero at high W
        NUMAX = (WMAX**2-M2+Q2)/TWOM
        LONUF2 = F2LONU(ITARG,Q2,NUMAX)
c        F2 = LONUF2 * EXP(10.0*(WMAX-W))
        F2 = LONUF2
      ENDIF

      IF (F2.LT.0.0) F2 = 0.0

      F2LOWW = F2

      RETURN
      END
       REAL FUNCTION F2LONU(ITARG,Q2EXT,NUEXT)
       IMPLICIT NONE
***************************************************************
* parametrizations for F2 at low nu    A V Kotwal 27dec93     *
* Original reference is F.W.Brasse et al, NP B39, 421 (1972)  *
* SLAC and DESY electroproduction and photoproduction         *
* data has been used in the fit. This fit can be used in the  *
* range 2<W<5 over a large Q2 range down to Q2=0              *
*                                                             *
* This parametrization is used in J. Franz et al,             *
* 'Inclusive Electron Scattering in the Low Q2 region'        *
* Z. Physics C 10 105-116 (1981)                              *
* in the range W>2, nu<6.5, 0.08<Q2<1.0                       *
*                                                             *
* ITARG=1 for H2, 2 for D2                                    *
* Q2EXT = Q2 in GeV^2                                         *
* NUEXT = Nu in GeV                                           *
***************************************************************
       REAL OMEGAW,NU,Q2,OMEGA,YOW,XPRIME,F2P,F2N,F2D,F2
       INTEGER ITARG
       REAL MW2,A2,B3,B4,B5,B6,B7
       PARAMETER (MW2=1.43,A2=0.42,B3=0.933,B4=-1.494,B5=9.021)
       PARAMETER (B6=-14.50,B7=6.453)
       REAL Q2EXT,NUEXT

       Q2 = Q2EXT
       NU = NUEXT
       IF (NU.LT.0.0) NU = 0.001
       IF (Q2.LT.0.0) Q2 = 0.001

       OMEGAW = (1.87654*NU+MW2)/(Q2+A2)
       OMEGA  =  1.87654*NU/Q2
       IF (OMEGA.LT.1.0) OMEGA = 1.0

       YOW    = 1.0 - 1.0/OMEGAW
       XPRIME = Q2/(1.87654*NU+0.8803)
       IF (XPRIME.GT.0.9999) XPRIME=0.9999

       F2P = B3+B4*YOW+B5*YOW**2+B6*YOW**3+B7*YOW**4
       F2P = OMEGAW*F2P*YOW**3/OMEGA
       IF (F2P.LT.0.0) F2P=0.0

       F2N = F2P*(1.0-XPRIME)
       IF (F2N.LT.0.0) F2N=0.0

       F2D = (F2N+F2P)/2.0

       IF (ITARG.EQ.1) THEN
         F2=F2P
       ELSEIF (ITARG.EQ.2) THEN
         F2=F2D
       ELSE
         F2=0.0
       ENDIF

       IF (F2.LE.0.0) F2 = 0.0

9999   CONTINUE
       F2LONU=F2
       RETURN
       END
       REAL FUNCTION F2PRES(Q2EXT,WEXT)
       IMPLICIT NONE
***********************************************************
* returns F2(proton) in the resonance region 1.1<W<2.0    *
* paremetrization of the cross-section is obtained from   *
* F.W.Brasse et al,                                       *
* ' Parametrization of the Q2 Dependence of the           *
* Gamma_v-P Total Cross-Sections in the Resonance Region' *
* Nuclear Physics B110, 413 (1976)                        *
*                                                         *
* for high energy muon scattering, the virtual photon     *
* polarization is close to 1 for low W, low Q2 scattering *
* such as resonance excitation.                           *
* Hence the parametrization for epsilon>0.9 is used from  *
* this reference. The average epsilon for this data is    *
* 0.957 . This is used together with R=sigmaL/sigmaT to   *
* convert the cross-section to F2                         *
*                                                         *
* A V Kotwal 28dec93                                      *
* inputs: Q2EXT = Q2 in GeV^2                             *
*         WEXT  = W  in GeV                               *
***********************************************************
       REAL M2,TWOM,F2P,Q2,NU,W2,W,Q,Q0,WLO,WHI,LNQOQ0,M,XX
       INTEGER NWBIN,IWBIN,N
       PARAMETER (NWBIN=56)
       REAL WMIN,WINTER,WMAX,DW1,DW2
       PARAMETER (WMIN=1.11,WINTER=1.77,WMAX=1.99)
       PARAMETER (DW1=0.015,DW2=0.02)
       PARAMETER (M=0.93828)    ! proton mass
       PARAMETER (M2=M*M,TWOM=2.0*M)
       REAL EPSILON,R,SIGMA,PI2AL4,Y,Y1,Y2,GD2,STRUW2,ERR
       PARAMETER (EPSILON=0.957)
       PARAMETER (PI2AL4=112.175)  ! 4.pi**2.alpha_em (microbarn.GeV**2)
       REAL Q2EXT,WEXT
       REAL A(NWBIN),B(NWBIN),C(NWBIN)
       DATA A /
     &  5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178
     & ,6.125,5.999,5.769,5.622,5.431,5.288,5.175,5.131
     & ,5.003,5.065,5.045,5.078,5.145,5.156,5.234,5.298
     & ,5.371,5.457,5.543,5.519,5.465,5.384,5.341,5.328
     & ,5.275,5.296,5.330,5.375,5.428,5.478,5.443,5.390
     & ,5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158
     & ,5.159,5.178,5.182,5.195,5.160,5.195,5.163,5.172 /
       DATA B /
     &  0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878
     & ,1.887,1.927,2.041,2.089,2.148,2.205,2.344,2.324
     & ,2.535,2.464,2.564,2.610,2.609,2.678,2.771,2.890
     & ,2.982,3.157,3.188,3.315,3.375,3.450,3.477,3.471
     & ,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519
     & ,4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322
     & ,5.546,5.623,5.775,5.894,6.138,6.151,6.301,6.542 /
       DATA C /
     &   0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080
     & ,-0.065,-0.056,-0.065,-0.056,-0.043,-0.034,-0.054,-0.018
     & ,-0.046,-0.015,-0.029,-0.048,-0.032,-0.046,-0.084,-0.115
     & ,-0.105,-0.159,-0.164,-0.181,-0.203,-0.220,-0.245,-0.264
     & ,-0.239,-0.302,-0.299,-0.318,-0.388,-0.393,-0.466,-0.588
     & ,-0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798
     & ,-1.048,-0.980,-1.021,-1.092,-1.313,-1.341,-1.266,-1.473 /

      Q2 = Q2EXT
      W  = WEXT
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (W.LT.M) W = M
      W2 = W*W
      NU = (W2+Q2-M2)/TWOM

      IF (NU.LT.0.001) NU = 0.001

* Q = three momentum transfer to the hadronic system in the lab frame
      Q = ABS(SQRT( Q2 + NU*NU ))
      IF (Q.LT.0.01) Q = 0.01

* Q0 is the value of Q at q2=0 for the same W
      Q0 = (W2 - M2)/TWOM
      IF (Q0.LT.0.01) Q0 = 0.01

c.....find the bin of W, and bin edges
      IF (W.GE.WMAX) THEN
        IWBIN = NWBIN
        WLO   = WMAX
        WHI   = WMAX
      ELSEIF (W.GT.WINTER) THEN
        N     = INT((W-WINTER)/DW2)
        IWBIN = N + 45
        WLO   = FLOAT(N)*DW2 + WINTER
        WHI   = WLO + DW2
      ELSEIF (W.GT.WMIN) THEN
        N     = INT((W-WMIN)/DW1)
        IWBIN = N + 1
        WLO   = FLOAT(N)*DW1 + WMIN
        WHI   = WLO + DW1
      ELSE
        IWBIN = 1
        WLO   = WMIN
        WHI   = WMIN
      ENDIF

      LNQOQ0 = ALOG(Q/Q0)

c.....eqn 3
c.....parameter d is fixed d=3 in section 3.1
      Y1 = A(IWBIN) + B(IWBIN)*LNQOQ0 + C(IWBIN)*(ABS(LNQOQ0))**3
      IF (Y1.LT.0.0) Y1 = 0.0

      IF (IWBIN.LT.NWBIN) THEN
       Y2=A(IWBIN+1)+B(IWBIN+1)*LNQOQ0+C(IWBIN+1)*(ABS(LNQOQ0))**3
c.....linear interpolation inside bin
       Y = Y1 + (W-WLO)*(Y2-Y1)/(WHI-WLO)
      ELSE
       Y = Y1
      ENDIF
      IF (Y.LT.0.0) Y = 0.0

* The parametrization returns log(sigma/GD**2) where GD is the
* nucleon dipole form factor, and sigma is the total virtual photoabsorption
* cross-section
      GD2   = 1.0 / (1.0 + Q2/0.71)**4
      SIGMA = GD2*EXP(Y)

c.....R=sigmaL/sigmaT
c.....use 1990 SLAC analysis result
      XX = Q2/(TWOM*NU)
      IF (XX.GT.0.99) XX=0.99
      IF (XX.LT.0.0) XX=0.0
      CALL R1990F(XX,Q2,R,ERR)

c.....'An introduction to Quarks and Leptons', F.E.Close
c......eqn 9.44,9.46
c......also eqn 1 from reference
      STRUW2 = SIGMA*(1.0+R)*Q2 / ((1.0+EPSILON*R)*(Q2+NU*NU))
c......the Hand convention for the virtual photon flux has been followed
c......(reference 3 in this reference)
      STRUW2 = STRUW2*Q0/PI2AL4
      F2P    = NU*STRUW2
      IF (F2P.LT.0.0) F2P = 0.0

9999  CONTINUE
      F2PRES = F2P

      RETURN
      END
            
              
!---------------------------------------------------------------------

       Subroutine F2SMC98(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the SMC parametrisation
! in SMC: Adeva et al (SMC) Phys. Rev. D 58, 112001 
! in Tables 12 and 13.


      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi,errd_lo,errd_hi,FNP_NMC,F2D
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2SMC98_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       CALL F2SMC98_DO(2,x4,qsq4,F2D,errd_lo,errd_hi)  !deuteron
       F2 = 2.*F2D- F2
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END

      SUBROUTINE F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)/
     >  -0.24997, 2.3963, 0.22896, 0.08498, 3.8608, -7.4143, 3.4342,  !p      
     >  -0.28151, 1.0115, 0.08415,-0.72973, 2.8647, -2.5328, 0.47477 /
      real*8 b(4,2)/ 0.11411, -2.2356, 0.03115, 0.02135,    !p
     >               0.20040, -2.5154, 0.02599, 0.01859/ 

      real*8 c(4,2)/ -1.4517, 8.4745, -34.379, 45.888,    !p
     >               -1.3569, 7.8938, -29.117, 37.657 /
              
!lower limits
      real*8 al(7,2)/
     >  -0.25196, 2.4297, 0.21913,  0.21630, 3.4645, -6.9887, 3.2771, !p  
     >  -0.28178, 1.1694, 0.09973, -0.85884, 3.4541, -3.3995, 0.86034/
      real*8 bl(4,2)/ 0.13074, -2.2465, 0.02995, 0.02039, 
     >                0.20865, -2.5475, 0.02429, 0.01760/  
      real*8 cl(4,2)/ -1.4715, 8.9108, -35.714, 47.338,
     >                -1.3513, 8.3602, -31.710, 41.106/

!upper limits
      real*8 au(7,2)/
     > -0.24810, 2.3632, 0.23643, -0.03241, 4.2268, -7.8120, 3.5822,    !p
     > -0.28047, 0.8217, 0.06904, -0.60191, 2.2618, -1.6507, 0.08909 /
      real*8 bu(4,2)/0.09734, -2.2254, 0.03239, 0.02233,
     >               0.18711, -2.4711, 0.02802, 0.01973  /
      real*8 cu(4,2)/ -1.4361, 8.1084, -33.306, 44.717,
     >                -1.3762, 7.6113, -27.267, 35.100  /

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER      
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
       Subroutine F2NMC(T,x,qsq,F2,err_f2) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in Phy Lett. b295 159-168 Proton and Deuteron F2 Structure Functions  
! in deep inelastic muon scattering   P. Amaudruz et. al.               
!                                                                       
                                                                        
                                                                        
      real a(7,2)/-0.1011,2.562,0.4121,-0.518,5.967,-10.197,4.685,      
     >            -0.0996,2.489,0.4684,-1.924,8.159,-10.893,4.535/      
      real b(4,2)/0.364,-2.764,0.0150,0.0186,                           
     >            0.252,-2.713,0.0254,0.0299/                           
      real c(4,2)/-1.179,8.24,-36.36,47.76,                             
     >            -1.221,7.50,-30.49,40.23/                             
      real Lam/0.25/                                                    
      real Qsqo/20.0/                                                   
      real log_term                                                     
      real A_x,B_x,C_x,F2,err_f2,x,qsq                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
                                                                        
      A_x=x**a(1,t)*(1-x)**a(2,t)*(a(3,t)+a(4,t)*(1-x)+a(5,t)*(1-x)**2  
     >    +a(6,t)*(1-x)**3+a(7,t)*(1-x)**4)                             
      B_x=b(1,t)+B(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = alog(qsq/lam**2)/alog(qsqo/lam**2)                     
      F2=A_x*(log_term)**B_x*(1+C_x/qsq)                                
      return                                                            
      end                                                               

!---------------------------------------------------------------------
      REAL FUNCTION FITEMC_N(X,A,Z,GOODFIT)                                         
!---------------------------------------------------------------------  
! Modified FITEMC.F with Neutron excess correction and proton=1 added 8/19/98
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! Z= number of protons
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 8/19/98 **  If proton, return 1.
! 8/19/98 **  Add neutron excess correction.
! 11/05 Modified PYB to return value at x=0.7 if x>0.7 since
!       Fermi smearing in INELASTU will account for Fermi smearing rise
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL ALPHA, C,LN_C,X,A,Z ,X_U,SIG_N_P,F_IS
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points
!Term    Coeficient     Error
      REAL*8 ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                 
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /                                  

             !
!Chisq=         22.    for 30 points
!Term    Coeficient     Error 
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        

      IF(A.LT.1.5) THEN    ! Added 8/19/98
       FITEMC_N=1.
       GOODFIT=.TRUE.
       RETURN
      ENDIF                                                                                
c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c      IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC_N  =  C *A**ALPHA    !isoscaler
      SIG_N_P = 1.-0.8*X_U
      F_IS = .5*(1.+ SIG_N_P)/(Z/A +(1.-Z/A)*SIG_N_P)
      FITEMC_N = FITEMC_N/F_IS
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o inelstu.o -c inelstu.f"; -*-
       SUBROUTINE INELSTU(X,q2,F2,W1,W2,IMOD)
! New call to Linda Stuarts fit to resonance Larry's DIS fit
! Steve Rock 10/14/93

*******************************************************************************
* This subroutine returns W1 and W2, the inelastic structure functions
*
* 2/93, LMS.
******************************************************************************
       IMPLICIT NONE

       INTEGER IMOD
       REAL X, Q2, NU, W1, W2,  MP/0.93828/, MP2/.8803/,
     >        WSQ, WW1, WW2, FRAC
       REAL   F2,R,DR,ST,SY,SLOPE,DSLOPE
       LOGICAL GOOD

       W1 = 0.0
       W2 = 0.0

       WSQ = MP2   + Q2*(1./X - 1. )
       NU = (WSQ +Q2 -MP2)/(2.*MP)
       IF(WSQ.LT.1.16) RETURN

! Get proton structure functions
       CALL R1998(X,Q2,R,DR,GOOD)
       IF(WSQ.GT.3.5) THEN
         CALL F2GLOB(X,Q2,'H',IMOD,F2,ST,SY,SLOPE,DSLOPE,GOOD)
         W2 = F2/NU
         W1 = (1.0 + NU*NU/Q2)*W2/(1.0 + R)
       ENDIF
! H2 model only officially works above q2=.75
       IF(Q2.LT.10.0.AND.WSQ.LT.4.3) THEN
         IF(WSQ.LE.3.5) THEN
           CALL H2MODEL(Q2,WSQ,W1,W2)
           F2 = NU*W2
         ELSE
           CALL H2MODEL(Q2,WSQ,WW1,WW2)
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + WW1*(1.0 - FRAC)
           W2 = W2*FRAC + WW2*(1.0 - FRAC)
         ENDIF
       ENDIF
       F2 = NU*W2
       RETURN
       END

        SUBROUTINE H2MODELLINDA(QSQ,WSQ,W1,W2)
!**** NO LONGER USED 5/16/06
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
********************************************************************************
        IMPLICIT NONE
        logical goroper/.false./

        REAL*4  WSQ, QSQ, W1, W2, R_NRES, SIG_RES, SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, DIPOLE, K, NU,
     >          TAU, PI/3.14159265/, AM/.93828/, ALPHA/0.00729735/,
     >          CONV/0.0025767/,sigroper(2)

! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.16) return
        IF(WSQ.GT.4.3.OR.QSQ.GT.10.0) THEN
          WRITE(6,'(''  H2MODEL called outside of kinematic range'')')
          RETURN
        ENDIF

! H2MOD_FIT returns transverse cross sections in units of
! microbarns/(dipole FF)**2
!        CALL H2MOD_FIT_old(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2)
        CALL H2MOD_FIT(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)

        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.

        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
        SIG_RES  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1)) *
     >    DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES
        SIGL = R_NRES*SIG_NRES(1)

! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END


      SUBROUTINE H2MOD_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J
      REAL W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, 
     >     RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER

      REAL PI/3.14159265/, ALPHA/0.00729735/, AM/.93828/,
     >     MPPI/1.07783/, MDELTA/1.2340/, MPI/0.13957/, 
     >     GAM1WID/0.0800/, GAM2WID/0.0900/, MASS1/1.5045/,
     >     ROPERWID/0.0500/, MASSROPER/1.4000/, 
     >     MASS2/1.6850/, DELWID/0.1200/, FIT(30), SQRTWDIF,
     >     XR/0.1800/, MQDEP/3.40/

      REAL RCOEF(25)/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >  -1.5416E+02,  -1.4891E+02,   2.4102E+02,   2.5823E+00, 
     >   7.1004E+00,  -8.9771E+00,   1.3744E+00,  -1.2085E+00, 
     >   1.1218E-01/

      REAL COEF(22)/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/

      REAL RERR1(200)/
     >  2.6120E+02,-9.4211E+02, 4.0844E+03, 7.4994E+02,-3.5317E+03,
     >  3.1703E+03, 2.1238E-01,-6.1568E-01, 4.1479E-01, 1.9720E-02,
     >  1.2891E-01,-4.1615E+00, 4.7246E+00, 2.8090E-03, 6.1657E-02,
     >  1.3120E+00,-9.4379E+00, 9.0902E+00,-1.3636E-03, 2.8054E-02,
     >  9.6123E-02,-1.1465E+03, 3.9099E+03,-3.0097E+03,-1.0604E+01,
     > -1.2214E+00,-8.3549E-01, 1.3696E+04, 3.9085E+03,-1.5369E+04,
     >  1.2745E+04, 2.9942E+01, 7.7268E+00, 1.0113E+01,-4.3868E+04,
     >  1.5709E+05,-3.0207E+03, 1.2809E+04,-1.1075E+04,-2.0442E+01,
     > -7.5843E+00,-1.0773E+01, 3.2597E+04,-1.2457E+05, 1.0283E+05,
     > -1.6960E-01, 5.9410E-01,-4.5486E-01,-1.0715E-02,-2.6512E-03,
     >  1.0153E-03, 6.4074E+00,-1.9189E+01, 1.3612E+01, 9.2486E-03,
     >  2.7904E-01, 6.3576E+00,-7.8552E+00, 1.5302E-02,-1.1506E-01,
     > -4.7552E-02,-1.0171E+01,-1.5884E+00, 1.6223E+01,-1.1379E-04,
     >  4.9212E-01,-2.4354E+00, 1.7921E+01,-1.7223E+01, 4.0778E-03,
     > -4.5558E-02,-1.8539E-01, 7.9930E+00,-7.1588E+01, 7.1512E+01,
     > -2.1529E-03, 1.8337E-01, 7.7590E-01, 7.3007E+02,-2.5219E+03,
     >  1.9547E+03, 6.1102E+00, 1.2970E+00,-1.3084E+00,-9.4932E+03,
     >  3.1542E+04,-2.3894E+04,-5.9583E+00, 8.1005E-02, 3.6885E-01,
     >  9.5708E+03,-2.4911E+03, 9.4342E+03,-7.7120E+03,-1.8608E+01,
     > -1.1065E+00, 6.5015E+00, 3.1755E+04,-1.1529E+05, 9.1964E+04,
     >  1.8347E+01,-2.5899E+00, 7.1169E-01,-3.2268E+04, 1.1891E+05,
     >  1.9339E+03,-7.7737E+03, 6.6128E+03, 1.3392E+01,-7.3587E-02,
     > -4.9353E+00,-2.4223E+04, 9.2768E+04,-7.6712E+04,-1.3210E+01,
     >  1.2513E+00,-4.5156E+00, 2.4541E+04,-9.5131E+04, 7.8848E+04,
     >  1.9952E-02,-7.1332E-02, 5.5522E-02, 9.8804E-04, 2.3682E-04,
     > -7.9762E-05,-6.3638E-01, 1.9492E+00,-1.4036E+00,-9.9312E-04,
     > -7.8634E-05, 8.2617E-05, 6.8002E-01,-2.1138E+00, 1.5308E+00,
     >  1.3008E-04,-1.0389E+02, 3.5942E+02,-2.7883E+02,-6.0671E-01,
     > -1.3016E-01, 1.4621E-01, 1.2841E+03,-4.3361E+03, 3.3132E+03,
     >  7.0701E-01, 1.2805E-01, 1.3355E-01,-1.4645E+03, 4.9522E+03,
     > -3.7686E+03,-1.0047E-01, 2.7406E+02, 3.5483E+02,-1.3433E+03,
     >  1.0978E+03, 1.9033E+00, 5.3726E-02,-8.1621E-01,-4.3612E+03,
     >  1.6110E+04,-1.2957E+04,-2.2247E+00,-2.1299E-01,-5.8178E-01,
     >  4.9755E+03,-1.8393E+04, 1.4724E+04, 3.1774E-01,-9.2555E+02,
     >  3.4086E+03,-2.7508E+02, 1.1025E+03,-9.3653E+02,-1.4100E+00,
     >  7.3163E-02, 6.6492E-01, 3.3590E+03,-1.3073E+04, 1.0893E+04,
     >  1.6311E+00, 2.4826E-01, 8.3308E-01,-3.7999E+03, 1.4772E+04,
     > -1.2252E+04,-2.3255E-01, 7.0167E+02,-2.7162E+03, 2.2434E+03,
     >  3.0688E+00,-1.0328E+01, 7.8828E+00, 3.6601E-03, 1.3367E-03,
     > -2.9672E-03,-3.2441E+01, 1.0979E+02,-8.3795E+01,-6.6345E-03/


      real rerr2(125)/
     > 3.7074E-02,
     >-5.7300E-02, 1.5212E-02, 4.5952E-04,
     > 1.1568E-04,-2.9315E-04,-4.6018E-01,
     > 9.3624E-01,-4.5908E-01,-6.2914E-05,
     > 1.1699E-03, 2.0141E-03, 6.9968E-02,
     >-1.9348E-01, 1.2176E-01, 5.4214E-07,
     > 1.3845E-04, 2.5311E-03,-2.5396E-03,
     >-1.2757E-04, 2.4971E-04,-1.2737E-04,
     > 7.2023E-03,-4.1404E-03, 4.6704E-04,
     > -4.6388E-03,-5.2545E-03, 4.0159E+01,-1.3481E+02, 1.0186E+02,
     >  1.1796E-03,-9.1088E+00, 3.0200E+01,-2.2552E+01, 4.3562E-01,
     > -1.0404E+01, 3.8414E+01,-3.0978E+01,-1.4730E-02, 4.6327E-03,
     >  1.9716E-02, 1.1236E+02,-4.1952E+02, 3.3862E+02, 2.4150E-02,
     >  1.1098E-02, 2.0122E-02,-1.3812E+02, 5.1058E+02,-4.0773E+02,
     > -4.1791E-03, 3.0702E+01,-1.1132E+02, 8.7622E+01,-1.4199E+00,
     >  5.0230E+00, 8.0171E+00,-3.1384E+01, 2.6350E+01, 1.3147E-02,
     > -6.1508E-03,-1.6808E-02,-8.7538E+01, 3.4530E+02,-2.8922E+02,
     > -1.9581E-02,-1.0895E-02,-2.4705E-02, 1.0611E+02,-4.1369E+02,
     >  3.4296E+02, 3.2847E-03,-2.3191E+01, 8.8502E+01,-7.2288E+01,
     >  1.0469E+00,-3.8873E+00, 3.1142E+00,
     > 1.1348E+00,-1.7657E+00, 4.7686E-01,
     > 1.6653E-02, 4.3488E-04,-7.5168E-03,
     >-1.6696E+01, 3.4692E+01,-1.7470E+01,
     >-4.9697E-03, 4.4232E-02, 5.7617E-02,
     > 5.7800E+00,-1.3886E+01, 7.9819E+00,
     > 3.4744E-04,-5.4411E-01, 1.2683E+00,
     >-7.0771E-01, 1.1282E-02,-2.4800E-02,
     > 1.2909E-02, 1.5171E-01,-6.0417E-01,
     > 7.7405E-01,-5.8981E-02,-5.8502E-03,
     > 8.8611E-04, 5.8326E-03, 6.5418E+00,
     >-1.2978E+01, 6.1069E+00, 1.2462E-03,
     >-1.8442E-02,-2.7954E-02,-1.8335E+00,
     > 4.3674E+00,-2.4393E+00,-6.2354E-05,
     > 1.4746E-01,-3.4127E-01, 1.8285E-01,
     >-3.0479E-03, 6.8138E-03,-3.4673E-03,
     >-7.5270E-02, 4.0914E-02/

      REAL ERR1(200)/
     >  3.7797E+02,-1.2732E+03, 4.8470E+03, 9.7589E+02,-3.9592E+03,
     >  3.3447E+03, 1.9629E-01,-4.2402E-01, 1.9757E-01, 3.0613E-02,
     > -4.0257E-01,-2.0922E+00, 3.0126E+00, 3.8385E-03, 7.3553E-02,
     >  1.4084E+00,-8.4718E+00, 7.8586E+00,-1.6484E-03, 2.2185E-02,
     >  7.4896E-02,-1.5627E+03, 5.0106E+03,-3.7125E+03,-1.1701E+01,
     > -6.9186E-01,-1.4263E+00, 1.5792E+04, 5.0288E+03,-1.7793E+04,
     >  1.3974E+04, 3.1643E+01, 5.0040E+00, 9.9958E+00,-4.8540E+04,
     >  1.6247E+05,-3.7498E+03, 1.4066E+04,-1.1461E+04,-2.0806E+01,
     > -5.0428E+00,-9.7813E+00, 3.5056E+04,-1.2382E+05, 9.7850E+04,
     > -2.0038E-01, 5.9769E-01,-4.0397E-01,-1.5776E-02,-3.7509E-03,
     >  5.7496E-04, 7.2218E+00,-2.0335E+01, 1.3722E+01, 1.2562E-02,
     >  1.4708E+00, 1.8510E+00,-4.1856E+00, 1.9572E-02,-1.3469E-01,
     > -3.7791E-02,-1.5215E+01, 1.8843E+01,-9.9384E-01, 5.4133E-04,
     >  5.6775E-01,-2.4158E+00, 1.5245E+01,-1.4180E+01, 5.3668E-03,
     > -3.5419E-02,-1.4360E-01, 7.8707E+00,-5.7677E+01, 5.5406E+01,
     > -7.5727E-04, 1.4127E-01, 5.8964E-01, 1.0277E+03,-3.3407E+03,
     >  2.4943E+03, 6.1372E+00, 2.0731E+00,-1.0628E-01,-1.1445E+04,
     >  3.6033E+04,-2.6376E+04,-6.4849E+00,-1.5437E+00,-3.1093E+00,
     >  1.1966E+04,-3.3062E+03, 1.1473E+04,-8.9323E+03,-1.7658E+01,
     > -3.0298E+00, 2.4862E+00, 3.6140E+04,-1.2237E+05, 9.3797E+04,
     >  1.8377E+01, 2.4649E-01, 9.5713E+00,-3.7362E+04, 1.2613E+05,
     >  2.4733E+03,-8.9836E+03, 7.2301E+03, 1.2133E+01, 1.0120E+00,
     > -2.0972E+00,-2.6581E+04, 9.4364E+04,-7.4804E+04,-1.2397E+01,
     >  5.8276E-01,-9.1893E+00, 2.7145E+04,-9.6250E+04, 7.6086E+04,
     >  2.4070E-02,-7.3772E-02, 5.1165E-02, 1.4597E-03, 3.3977E-04,
     > -2.6275E-05,-7.2542E-01, 2.0676E+00,-1.4052E+00,-1.3577E-03,
     > -1.4477E-04,-8.5451E-05, 7.4811E-01,-2.1217E+00, 1.4288E+00,
     >  1.7439E-04,-1.6022E+02, 5.2231E+02,-3.9172E+02,-4.1771E-01,
     > -2.3133E-01,-1.9119E-02, 1.6931E+03,-5.4146E+03, 4.0099E+03,
     >  6.5228E-01, 4.5766E-01, 6.7254E-01,-2.0266E+03, 6.3551E+03,
     > -4.6404E+03,-9.4689E-02, 4.2768E+02, 5.1531E+02,-1.7829E+03,
     >  1.3890E+03, 1.1798E+00, 3.1335E-01,-2.5902E-01,-5.3955E+03,
     >  1.8502E+04,-1.4311E+04,-1.8045E+00,-9.6753E-01,-2.0260E+00,
     >  6.3626E+03,-2.1445E+04, 1.6387E+04, 2.6350E-01,-1.3456E+03,
     >  4.5055E+03,-3.8598E+02, 1.3911E+03,-1.1170E+03,-7.9328E-01,
     > -7.6840E-02, 2.5967E-01, 4.0005E+03,-1.4347E+04, 1.1455E+04,
     >  1.1795E+00, 6.2629E-01, 1.6961E+00,-4.6485E+03, 1.6399E+04,
     > -1.2954E+04,-1.7187E-01, 9.8638E+02,-3.4363E+03, 2.7002E+03,
     >  6.0266E+00,-1.9528E+01, 1.4686E+01,-1.7956E-02, 3.3364E-03,
     >  1.2080E-03,-5.5018E+01, 1.7933E+02,-1.3517E+02, 7.9955E-03/


      REAL ERR2(53)/
     > -2.1546E-02,-2.3493E-02, 7.4315E+01,-2.3518E+02, 1.7398E+02,
     > -6.4429E-04,-1.9950E+01, 6.3147E+01,-4.6881E+01, 1.2816E+00,
     > -1.9366E+01, 6.5755E+01,-5.0971E+01, 5.7005E-02, 3.3439E-04,
     >  5.5786E-03, 1.7715E+02,-6.1369E+02, 4.7999E+02,-2.9558E-02,
     >  5.5461E-02, 7.1075E-02,-2.3560E+02, 7.9154E+02,-6.0792E+02,
     >  2.7242E-03, 6.3265E+01,-2.0981E+02, 1.6050E+02,-4.0749E+00,
     >  1.3388E+01, 1.4562E+01,-5.1058E+01, 4.0557E+01,-4.3474E-02,
     > -4.4868E-03,-6.3041E-03,-1.3274E+02, 4.7814E+02,-3.8441E+02,
     >  2.5677E-02,-3.8538E-02,-5.8204E-02, 1.7424E+02,-6.0799E+02,
     >  4.8014E+02,-2.6425E-03,-4.6992E+01, 1.6058E+02,-1.2570E+02,
     >  3.0554E+00,-1.0258E+01, 7.9929E+00/

      LOGICAL FIRST/.TRUE./

! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
*       write(*,*)'dd1',mdelta,mpi,epircm,first
!Define error matrix:
        K = 0
        if (goroper) then 
          DO J = 1,25
            DO I = 1,J
              K = K + 1
              if (k.le.200) RERRMAT(I,J) = RERR1(K)
              if (k.le.325.and.k.gt.200) RERRMAT(I,J)=RERR2(K-200)
            ENDDO
          ENDDO
        endif
       if (.not.goroper) then  
          DO J = 1,22
            DO I = 1,J
              K = K + 1
              if (k.le.200) ERRMAT(I,J) = ERR1(K)
              if (k.le.253.and.k.gt.200) ERRMAT(I,J)=ERR2(K-200)
            ENDDO
          ENDDO
        endif
        if (goroper) then
          DO J = 1,25
              DO I = J+1,25
                RERRMAT(I,J) = RERRMAT(J,I)
              ENDDO
          ENDDO
        endif
        if (.not.goroper) then
          DO J = 1,22
              DO I = J+1,22
                ERRMAT(I,J) = ERRMAT(J,I)
              ENDDO
          ENDDO
        endif
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN
*      write(*,*)'sigdel=',sigdel,mdelta

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(23) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
      ENDDO
      if (goroper) FIT(24)  = FIT(23)/sqrt(Q2)
      if (goroper) FIT(25)  = FIT(23)/q2
      DO I = 1,4
        FIT(I + 12)  = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 16)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 19)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,25
          RSIGERR(J) = 0.0
          DO I = 1,25
            RSIGERR(J) = RSIGERR(J) + FIT(J)*FIT(I)*RERRMAT(I,J)
            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*RERRMAT(I,J)
          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
             SIG_RES2(2) = SIG_RES2(2) + RSIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
             SIG_RES1(2) = SIG_RES1(2) + RSIGERR(J)
          elseIF(j.ge.23.and.j.le.25) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
            SIGroper(2) = SIGroper(2) + RSIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
            SIG_NRES(2) = SIG_NRES(2) + RSIGERR(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
          SIGERR(J) = 0.0
          DO I = 1,22
            SIGERR(J) = SIGERR(J) + FIT(J)*FIT(I)*ERRMAT(I,J)
            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*ERRMAT(I,J)
          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
             SIG_RES2(2) = SIG_RES2(2) + SIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
            SIG_RES1(2) = SIG_RES1(2) + SIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
            SIG_NRES(2) = SIG_NRES(2) + SIGERR(J)
          ENDIF
        ENDDO
      endif

! ERRCHECK should agree with SIGTOTER.
C      ERRCHECK = SQRT(ABS(SIG_RES2(2) + SIG_RES1(2) + SIG_NRES(2)))
C      SIGTOTER = SQRT(SIGTOTER)
      SIG_RES2(2) = SQRT(ABS(SIG_RES2(2)))
      SIG_RES1(2) = SQRT(ABS(SIG_RES1(2)))
      SIG_NRES(2) = SQRT(ABS(SIG_NRES(2)))

      RETURN
      END

      SUBROUTINE H2MOD_FIT_old(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2)
********************************************************************************
* This is a 22 parameter model fit to the NE11 and E133 data and three
* QSQ points of generated BRASSE data to constrain fits at low QSQ. It has
* three background terms and three resonances. Each of these is multiplied
* by a three-term polynomial in Q2. The DELTA resonance has a four term
* polynomial.
*
* 8/91, LMS.
********************************************************************************
      IMPLICIT NONE

      INTEGER I, J
      REAL W2, Q2, SIG_NRES, SIG_RES1, SIG_RES2,  KLR,
     >     KRCM/0./,EPIRCM/0./, PPIRCM, GG, GPI, DEN, SIGDEL,
     >     W, KCM, K, EPICM, PPICM, WDIF

      REAL PI/3.14159265/, ALPHA/0.00729735/, AM/.93828/,
     >     MPPI/1.07783/, MDELTA/1.2355/, MPI/0.13957/,
     >     GAM1WID/0.0676/, GAM2WID/0.1025/, MASS1/1.5062/,
     >     MASS2/1.6800/, DELWID/0.1260/, FIT(30), SQRTWDIF,
     >     XR/0.2060/, MQDEP/2.07/

      REAL COEF(22)/
     >  -0.3982E+04,   0.1764E+05,  -0.1590E+05,   0.1391E+02,
     >  -0.1233E+02,  -0.8206E+01,   0.9006E+04,  -0.3581E+05,
     >   0.3094E+05,  -0.2668E+01,   0.3306E+02,   0.2652E+02,
     >  -0.3134E+04,   0.2059E+05,  -0.1817E+05,   0.1327E+00,
     >   0.6619E+03,  -0.3914E+04,   0.3579E+04,  -0.4032E+02,
     >   0.2233E+03,  -0.2066E+03/

      LOGICAL FIRST/.TRUE./

! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
      ENDIF
      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2)
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model.
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 +
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
      ENDDO
      DO I = 1,4
        FIT(I + 12)  = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 16)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 19)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO

! Find sig_t (in microbarns/gd**2).
      SIG_NRES = 0.0
      SIG_RES1 = 0.0
      SIG_RES2 = 0.0
      DO J = 1,22
        IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
          SIG_RES2 = SIG_RES2 + FIT(J)*COEF(J)
        ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
          SIG_RES1 = SIG_RES1 + FIT(J)*COEF(J)
        ELSE
          SIG_NRES = SIG_NRES + FIT(J)*COEF(J)
        ENDIF
      ENDDO

      RETURN
      END
!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
       Subroutine F2NMC_new(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in CERN_PPE/95-138  Sept 4, 1995  Proton and Deuteron F2 Structure Functions 
! in deep inelastic muon scattering   P. Arneodo et. al.               
!  Published in Phys.Lett.B364:107-115,1995 
!   e-Print Archive: hep-ph/9509406 
!                   
      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi ,FNP_NMC
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2NMC_NEW_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       F2 = F2 * FNP_NMC(X4,QSQ4)
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END

      SUBROUTINE F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)
     >  /-0.02778,  2.926,   1.0362, -1.840, 8.123, -13.074, 6.215, 
     >   -0.04858,  2.863,    .8367, -2.532, 9.145, -12.504, 5.473/ 
      real*8 b(4,2)/ 0.285,   -2.694,   0.0188,  0.0274,  
     >              -0.008,   -2.227,   0.0551,  0.0570/
      real*8 c(4,2)/-1.413,    9.366, -37.79,   47.10, 
     >              -1.509,    8.553, -31.20,   39.98/

!lower limits
      real*8 al(7,2)
     >  /-0.01705,  2.851,   0.8213, -1.156, 6.836, -11.681, 5.645, 
     >   -0.02732,  2.676,    .3966, -0.608, 4.946,  -7.994, 3.686/ 
      real*8 bl(4,2)/ 0.325,   -2.767,   0.0148,  0.0226,  
     >                0.141,   -2.464,   0.0299,  0.0396/
      real*8 cl(4,2)/-1.542,   10.549, -40.81,   49.12, 
     >               -2.128,   14.378, -47.76,   53.63/

!upper limits
      real*8 au(7,2)
     >  /-0.05711,  2.887,   0.9980, -1.758, 7.890, -12.696, 5.992, 
     >    -0.04715,  2.814,    .7286, -2.151, 8.662, -12.258, 5.452/ 
      real*8 bu(4,2)/ 0.247,   -2.611,   0.0243,  0.0307,  
     >               -0.048,  -2.114,   0.0672,  0.0677/
      real*8 cu(4,2)/-1.348,    8.548, -35.01,   44.43, 
     >               -1.517,    9.515, -34.94,   44.42/
      

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               

!---------------------------------------------------------------------
      REAL FUNCTION FNP_NMC(X,QSQ)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL X,QSQ,A,B,X2,X3
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      FNP_NMC = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END

! correction routine from Blok/Vlados for 99118 use only
	SUBROUTINE F2corr1H(x,fact_cor)
***** theta in degree
        implicit none
      COMMON /KIN/   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2

 	REAL DIF1,DIF2,DIF3,DIF4,x,fact_cor,th1,th2,diffac
	REAL*4 F1,F2,DIFYGL,DIFTHETA
	REAL*4 GIPOT,COSYGL,GIPOSMALL
	integer i

	real P1(2)
	real P2(2)
	real P3(2)
	real P4(2)
 
*============ "START" THETA=10.60 ========================
 
 	 p1(1)=0.97094 
	 p1(2)=-0.98026E-01 
         dif1=p1(1)+p1(2)*x 
         	
*============ "FINISH" THETA=10.60 ========================	
	
*============ "START" THETA=14.60 ========================
		 
	 p2(1)=0.99776 
	 p2(2)=-0.34140E-01 
         dif2=p2(1)+p2(2)*x     
         	
*============ "FINISH" THETA=14.60 ========================	
	
*============ "START" THETA=18.60 =========================

 	 p3(1)=0.99493 
	 p3(2)=0.30907E-01 
         dif3=p3(1)+p3(2)*x  
         
*============ "FINISH" THETA=18.60 ========================	
	
*============ "START" THETA=22.60 =========================

 	 p4(1)=1.0284 
	 p4(2)=-0.47091E-01 
         dif4=p4(1)+p4(2)*x     
     
      
*============ "FINISH" THETA=22.60 ========================

 
         if(theta.lt.14.60) then
	      f1=dif1
	      f2=dif2
	      th1=10.60
	      th2=14.60
	 endif
	 
         if(theta.ge.14.60.and.theta.lt.18.60) then
	      f1=dif2
	      f2=dif3
	      th1=14.60
	      th2=18.60
	 endif

         if(theta.ge.18.60) then
	      f1=dif3
	      f2=dif4
	      th1=18.60
	      th2=22.60
	 endif
	       
	 if(f1.le.f2) difygl=abs(theta-th1)
	 if(f1.gt.f2) difygl=abs(theta-th2)

	 diffac=abs(f1-f2)	 
	 diftheta=th2-th1	 
         gipot=sqrt(diftheta**2+diffac**2)
	 cosygl=diftheta/gipot	 	   
         giposmall=difygl/cosygl
	 fact_cor=sqrt(giposmall**2-difygl**2)	 

	 if(f1.lt.f2) THEN
	  fact_cor=fact_cor+f1
	   else
	  fact_cor=fact_cor+f2
	 ENDIF   
	  
	 	          
! 	 WRITE(*,*)'x=:Factor=',x,fact_cor
  	 
	 do i=1,2
	  P1(i)=0.
	  P2(i)=0.
	  P3(i)=0.
	  P4(i)=0.
	 enddo
	 
	 dif1=0.
	 dif2=0.
	 dif3=0.
	 dif4=0.
	 f1=0.
	 f2=0.
	 
	 END


	SUBROUTINE F2corr2H(x,fact_cor)
***** theta in degree
        implicit none

 	REAL DIF1,DIF2,DIF3,DIF4,x,fact_cor,th1,th2,diffac
	REAL*4 F1,F2,DIFYGL,DIFTHETA
	REAL*4 GIPOT,COSYGL,GIPOSMALL
      COMMON /KIN/   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,THeta,THR,Xx,Q2,EPS,W2,Y,ANU,SIN2
        integer i

	real P1(2)
	real P2(2)
	real P3(2)
	real P4(2)
 
*============ "START" THETA=10.60 ========================
      	 
	 p1(1)=0.94718 
	 p1(2)=-0.15652E-01 
         dif1=p1(1)+p1(2)*x 
          	
*============ "FINISH" THETA=10.60 ========================	
	
*============ "START" THETA=14.60 ========================
	 
	 p2(1)=1.0024  
	 p2(2)=-0.30468E-01  
         dif2=p2(1)+p2(2)*x   
   	
*============ "FINISH" THETA=14.60 ========================	
	
*============ "START" THETA=18.60 =========================

	 
	 p3(1)=1.0327 
	 p3(2)=-0.32571E-01 
         dif3=p3(1)+p3(2)*x   
      
       
*============ "FINISH" THETA=18.60 ========================	
	
*============ "START" THETA=22.60 =========================

 	 p4(1)=1.0846 
	 p4(2)=-0.13985 
         dif4=p4(1)+p4(2)*x   
  
*============ "FINISH" THETA=22.60 ========================

 
         if(theta.lt.14.60) then
	      f1=dif1
	      f2=dif2
	      th1=10.60
	      th2=14.60
	 endif
	 
         if(theta.ge.14.60.and.theta.lt.18.60) then
	      f1=dif2
	      f2=dif3
	      th1=14.60
	      th2=18.60
	 endif

         if(theta.ge.18.60) then
	      f1=dif3
	      f2=dif4
	      th1=18.60
	      th2=22.60
	 endif
	       
	 if(f1.le.f2) difygl=abs(theta-th1)
	 if(f1.gt.f2) difygl=abs(theta-th2)

	 diffac=abs(f1-f2)	 
	 diftheta=th2-th1	 
         gipot=sqrt(diftheta**2+diffac**2)
	 cosygl=diftheta/gipot	 	   
         giposmall=difygl/cosygl
	 fact_cor=sqrt(giposmall**2-difygl**2)	 

	 if(f1.lt.f2) THEN
	  fact_cor=fact_cor+f1
	   else
	  fact_cor=fact_cor+f2
	 ENDIF   
			  
!* 	 WRITE(*,*)'x=:Factor=',x,fact_cor
  	 
	 do i=1,2
	  P1(i)=0.
	  P2(i)=0.
	  P3(i)=0.
	  P4(i)=0.
	 enddo
	 
	 dif1=0.
	 dif2=0.
	 dif3=0.
	 dif4=0
	 f1=0.
	 f2=0. 
	    	 
	 END

       SUBROUTINE D2MODEL_IOANA(QSQ,WSQ,W1,W2)
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
*****************************************************************************
        IMPLICIT NONE

        REAL*4  WSQ, QSQ
        REAL*4  R_NRES, SIG_RES(2), SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, 
     >          DIPOLE, K, NU,TAU, PI, AM, ALPHA, W1, 
     >          W2,CONV,sig_roper(2)
	integer i
	real*4  xval(37)
	logical goroper/.false./
*
        COMMON/ioana/xval
*       
*
* N.I. ...
*

        DATA PI/3.14159265/, AM/.93828/, ALPHA/0.00729735/,
     >          CONV/0.0025767/

*
! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.17) return
        IF(WSQ.LT.1.17.OR.WSQ.GT.5.OR.
     >    QSQ.GT.10.0) THEN
           do i=1,2	
              SIG_NRES(i) = 0.0
              SIG_RES1(i) = 0.0
              SIG_RES2(i) = 0.0
              SIG_ROPER(i)= 0.0 
           enddo
!          WRITE(*,*)'H2MODEL_IOANA called outside of kinematic range'
          RETURN
        ENDIF

!  returns transverse cross sections in units of
! microbarns/(dipole FF)**2
        CALL i_d2_model(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     &                  SIG_ROPER,goroper,xval)
*        write(*,*)'1',SIG_NRES
        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.
        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
*        write(*,*)'1',SIG_NRES(1),DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1)+sig_roper(1))*
     &                 DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES(1)
        SIGL = R_NRES*SIG_NRES(1)
*	write(*,*) 'I am here'

!        write(33,*)SIG_NRES(1),SIG_RES(1),R_NRES,DIPOLE,sigt
! 33     format(F15.3,2x,F15.3,2x,F15.3,2x,F15.3,2x,F15.3)      

! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)
*        write(*,*)'from program',w1,w2 

!*** PYB DIVIDE BY 2 TO GET PER NUCLEON
        W1=W1/2.
        W2=W2/2.

        RETURN
        END
         SUBROUTINE i_d2_model(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper,XVAL)
********************************************************************************
*
* This subroutine calculates model cross sections for deuterium in resonance region. 
* Cross section is returned in nanobarns. 
* For all the resonance regions we use nonrelativistc
* Breit-Wigner shapes (L.Stuart). 
* This subroutine is based on h2model.f from SLAC
* 
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
********************************************************************************
      IMPLICIT NONE
*
      REAL*4 XVAL(37)
* 
      logical goroper
      INTEGER I, J
*
* modified by I.N. 04/09/97
*

       
      REAL*4 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     KRCM, EPIRCM, PPIRCM,   
     >     W, KCM, K, EPICM, PPICM, WDIF 
 
      REAL*4 PI, ALPHA, AM,
     >     MPPI, MPI,  
     >     XR  

      integer jj
      real*4 kr,meta
      real*4 mrho,mn,am2
      real*4 gamma_gamma,gamma_pi,denom_pi
      real*4 sigma_delta,sigma_second,sigma_third
*
      real*4 mass(3),width(3),qdep(3)
      real*4 coeff(3,0:3),nres_coeff(0:4,3)
* Define vectors : MASS(3), WIDTH(3), which will contain the masses and widths of
* the three resonances: 1=DELTA, 2=second res reg, 3=third res reg
*
      DATA PI/3.14159265/
      DATA ALPHA/0.00729735/
      DATA AM/.93828/
      DATA MPI/0.13957/
      DATA MPPI/1.07783/ 
      data Meta/0.54745/ 
      data Mrho/0.7681/
*      DATA XR/0.18/
	data Mn/0.938/
      Am2 = Am*Am
*
* The fit parameters are called  XVAL 
* Assign them to masses, widths, etc.
*
      xr=XVAL(1)
      do i=1,3
         mass(i)  = XVAL(i+1)
*         write(*,*)xval(i)
      enddo
      do i=1,3
         width(i)  = XVAL(i+4)
*         write(*,*)xval(i+3)
      enddo
      do i=1,3
         qdep(i)   = XVAL(i+7)
*         write(*,*)xval(i+6)
      enddo     
      k=11
      do i=1,3
         do j=0,3
            coeff(i,j) = XVAL(k)
*            write(*,*)xval(k)
            k=k+1
         enddo
      enddo
      do i=0,4
         do j=1,3
           nres_coeff(i,j) = XVAL(k)
***!!!            write(*,*)'3',xval(k)
           k=k+1
         enddo
      enddo
*      write(*,*)'Done assigning parameters'
*
*      write(*,*) 'w2=', w2
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
****************************************************************
* This part is not used for deuterium.
****************************************************************
* K equivalent real photon energy needed to excite resonance W
* K is in Lab frame, KCM is in CM frame
* 
         K = (W2 - Am2)/2./Am
         KCM = (W2 - Am2)/2./W
         KR = (Mass(1)*Mass(1)-Am2)/2./Am
         KRCM = (Mass(1)*Mass(1) - Am2)/2./Mass(1)
* Decay mode for Delta: N,pi:
********
      EPICM = 0.5*( W2 + Mpi*Mpi - Am2 )/W

      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      EPIRCM = 0.5*(Mass(1)*Mass(1) + Mpi*Mpi - Am2 )/Mass(1)   
      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      
********
* Now calculate the partial widths:
*  gamma_pi(i)
*      write(*,*)'mass(1),w2 ',mass(1),' ' ,w2
*      write(*,*)'width(1) ',width(1)
*      write(*,*)'PPICM, PPIRCM ',PPICM,' ',PPIRCM
*      write(*,*)'XR ' ,xr
      gamma_pi = width(1)*(PPICM/PPIRCM)**3*(PPIRCM*PPIRCM+XR*XR)/
     >          (PPICM*PPICM+XR*XR)
*
*      write(*,*)'iuhu 2'
      gamma_gamma = width(1)*(KCM/KRCM)**2*(KRCM*KRCM+XR*XR)/
     >          (KCM*KCM+XR*XR)

      denom_pi = (W2 - Mass(1)*Mass(1))**2 + (Mass(1)*gamma_pi)**2
*      write(*,*)'iuhu 3'
***********************************************************************
***********************************************************************
* Now calculate cross sections for Delta:
      sigma_delta = width(1)/((W-Mass(1)*(1. + Q2*qdep(1)))**2 
     >               + 0.25*width(1)*width(1))
* For the second and third resonance regions:
      sigma_second = width(2)/((W-Mass(2)*(1. + Q2*qdep(2)))**2 
     >               + 0.25*width(2)*width(2))

*      write(*,*)'width(2)',' ',width(2)
*      write(*,*)'mass(2)',' ',mass(2)

*      write(*,*)'sigma_second',' ',sigma_second

      sigma_third  = width(3)/((W-Mass(3)*(1. + Q2*qdep(3)))**2
     >               + 0.25*width(3)*width(3))

*      write(*,*)'sigma_third',' ',sigma_third

*      write(*,*)'width(3)',' ',width(3)
*      write(*,*)'mass(3)',' ',mass(3)

*      write(*,*)'iuhu 5'
*
* Put in the Q2 dependence of AH(Q2):
* AH(Q2) = coeff(i,j)*Q2**j
*
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(1) = 0.0
      do j=0,3  
         sig_res1(1) = sig_res1(1)+ sigma_delta*coeff(1,j)*Q2**j
         sig_res2(1) = sig_res2(1)+
     >               sigma_second*coeff(2,j)*Q2**j+
     >               sigma_third *coeff(3,j)*Q2**j
      enddo   
*      write(*,*)'sigma_d,sigma_23',' ',sig_res1(1),' ',sig_res2(1)
      do j=0,4
         do jj=1,3
            sig_nres(1) = sig_nres(1)+
     >          nres_coeff(j,jj)*q2**j*sqrt(wdif**(2*jj-1))
*          write(*,*)'2',sig_nres(1),nres_coeff(j,jj),j,wdif,jj
         enddo
      enddo     
*
      RETURN
      END

      REAL FUNCTION  ROCK_RES9(TARG,X,Q2)
! April 22 2003 Version of deuterium fit.
! No JLAB data

      IMPLICIT NONE
      INTEGER TARG ! 1=H, 2=D
      REAL Q2
      REAL X,F2INEL
      REAL*8 WM,B03
      CHARACTER*1 TARG_STR ! D or H      
      REAL MP2/.8803/,MP/.93828/ 
! DATE =2003/ 4/22  TIME=13: 3
! MODEL# 9
! HYDROGEN
      REAL*8 CH(40)/
     >  1.074100000, 0.279494473, 4.712278073, 1.739263857, 1.431155682,
     >  0.184080601, 1.225717385, 0.110613494, 0.214285318, 1.510392164,
     >  0.105595846, 0.309575075, 1.725500000, 0.138478762,-1.452323856,
     >  1.953800000, 0.342265136,-0.044517487, 0.093869839,-0.029404837,
     >  1.199536964, 0.671722156, 2.033898985, 0.024107348, 0.452992650,
     >  0.624721200, 0.054377925, 0.099384033, 0.897445331,13.205228916,
     >  6.284002261, 0.013794386, 0.030764850,-0.007500832, 0.025566716,
     >  0.016759263, 0.966305525, 1.058950737, 0.000000000, 0.000000000/
     > 
!deuterium
      REAL*8 CD(40)/
     >   1.05220000,  0.75530000,  3.30000000,  1.70000000,  3.50000000,
     >   1.04000000,  1.22000000,  0.10000000,  0.48000000,  1.50000000,
     >   0.08000000,  0.65000000,  1.70000000,  0.12000000,  0.74000000,
     >   1.90000000,  0.19000000, -0.17000000,  0.00960000,  0.03500000,
     >   3.50000000, -0.60000000,  4.70000000,  0.41100000,  0.41100000,
     >   0.10000000,  0.41100000,  0.10000000,  0.10000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000/
    
      IF(TARG.EQ.1) TARG_STR ='H'
      IF(TARG.EQ.2) TARG_STR ='D'
      WM = SQRT(MP2 +Q2*(1./X -1.) )

      CALL F2GLOB_FAST(X,Q2,TARG_STR,9,F2INEL)

      IF(TARG.EQ.1) THEN  ! H2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CH)   
      ELSE     ! D2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CD)   
      ENDIF
      RETURN
      END

!----------------------------------------------------------------

      REAL*8 FUNCTION B03(WM,QSQ,NC2,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(80),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,BRES
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.8803/, PM2/1.876512/, PM/.93828/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      INTEGER   NC2 !offset for coeficients (0 for h2)
      DATA      LSPIN/1,2,3,2/                                       

C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(NC2+22)+C(NC2+23)*(X-C(NC2+24))**2                                 
      PIEMSQ = (C(NC2+1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(NC2+1)) B1 = C(NC2+2)                                         
      EB1 = C(NC2+3)*(WM-C(NC2+1))                                              
      IF( (EB1.LE.25.).AND.(B1.NE.0.)) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(NC2+4)) B2 = (1.-C(NC2+2))                                    
      EB2 = C(NC2+5)*(WSQ-C(NC2+4)**2)                                          
      IF( (EB2.LE.25.0).AND.(B2.NE.0.)) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(NC2+2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(NC2+INDEX)                                            
           IF (I.EQ.1) RAM=C(NC2+INDEX)+C(NC2+18)*SQRT(QSQ)! +C(NC2+30)*QSQ**2
     >       + C(NC2+25)/(QSQ+C(NC2+26))
           RMA    = C(NC2+INDEX+1)
           IF (I.EQ.2) RAM =RAM + C(NC2+27)/(C(NC2+28)+QSQ)
           IF (I.EQ.3) THEN
              RMA=RMA*(1.D0 +C(NC2+20)/(1.D0+C(NC2+21)*QSQ))
              RAM =RAM + C(NC2+19)/(C(NC2+29)+QSQ) 
           ENDIF
           IF (I.EQ.4) RAM =RAM + C(NC2+30)/(C(NC2+31)+QSQ)
           RWD    = C(NC2+INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
                                                                        
C FORM VW2/F2                                                           
                                                                        
      B03 = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
                                                                        
      RETURN                                                            
      END    

      subroutine christy(W2,Q2,F1,R)
!------------------------------------------------------------------------
! subroutine to return proton structure function F1 and ratio R=sigl/sigt
!
! inputs are electron missing mass squared W2 (GeV**2) 
!            momentum trasfer squared Q2 (GeV**2)
! inputs and outputs are Real*8
! the file christy.dat is needed to use this subroutine
! Note: if W2<1.155 GeV**2, values of zero are returned (below threshold)
! Fit done by Eric Christy 11/04
! Reference this fit as ???       
!------------------------------------------------------------------------
      IMPLICIT NONE

      real*8 w2,q2,xval1(41),xvall(41),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F1,FL,R
      integer i
      logical first/.true./
 
      mp = .93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.

      if(first) then
        first=.false.
        open(unit=15,file='christy.dat',status='old')
        do i=1,41
          read(15,*) temp(1)  ! par #
          read(15,*) XVAL1(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        do i=1,41
          read(15,*) temp(1)  ! par #
          read(15,*) XVALL(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(15)
      endif

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2/(w2+q2-mp2)
      if(xb.le.0.0) return

      call resmod(1,w2,q2,XVAL1,F1)
      call resmod(2,w2,q2,XVALL,FL)
      if(F1.le.0.0) return

       if(F1.le.0.0) return

      F1 = F1/8.d0/pi/pi/alpha/1000.d0
      F1 = F1*(w2-mp2)
      FL = FL/8.d0/pi/pi/alpha/1000.d0
      FL = FL*(w2-mp2)*2.d0*xb
      R = FL/ (2.D0 * XB * F1)

      return
      end

      SUBROUTINE RESMOD(sf,w2,q2,xval,fn) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,fn,xval(41),mass(4),width(4)
      REAL*8 height(4),fn_del,fn_s11,fn_f15,rescoef(4,3)
      REAL*8 nr_coef(4,4),wdif,fn_nr,fn_4,w2temp,wtemp
      REAL*8 roper_mass,roper_width,roper_height
      REAL*8 roper_mparm,roper_exp,mq2(4)
      REAL*8 alpha,pi
      INTEGER i,j,num,sf


      mp = 0.93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036
      W = sqrt(w2)
      w2temp = w2 - .93828*.93828
      wtemp = sqrt(w2temp)
      wdif = w - (0.937272 + 0.137)
      xb = q2/(q2+w2-mp2)
      if(sf.EQ.1) xb = 1.d0
      fn_nr = 0.0d0
      do i=1,4
        height(i) = 0.
      enddo
      
      do i=1,4
        mass(i) = xval(i)
      enddo

      do i=1,4
        mq2(i) = xval(36 + i)
      enddo

      if(q2.GT.0.1) then

        mass(2) = mass(2)*exp(-(q2-0.1)/mq2(2)) 
     &            + mq2(1)*(1.-exp(-(q2-0.1)/mq2(2))) 

        mass(3) = mass(3)*exp(-(q2-0.1)/mq2(4)) 
     &            + mq2(3)*(1.-exp(-(q2-0.1)/mq2(4)))

      endif   
      do i=1,4
        width(i) = xval(4+i)
      enddo
      num = 0
      do i=1,4
       do j=1,3
         num = num + 1
         rescoef(i,j)=xval(8 + num)
c           write(6,*) i,j,num,rescoef(i,j)
c           height(i) = height(i)+rescoef(i,j)*q2**(-1.*(float(j-1)))
         enddo
         height(i) = rescoef(i,1)*
     &             (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
c         if(w2.LT.1.35) height(i) = height(i)/(w2-1.30)    
      enddo
      num = 0     
      do i=1,4
       do j=1,4
         num = num + 1
         nr_coef(i,j)=xval(20 + num)
c         write(6,*) i,j,num,nr_coef(i,j)         
       enddo
      enddo

      do i=1,5
       roper_mass = xval(37)
       roper_width = xval(38)
       roper_height = xval(39)
       roper_mparm = xval(40)
       roper_exp = xval(41)
      enddo
c      write(6,*) "constant coef's are:  ",height

CC   Calculate Breit-Wigners for the 3 resonance regions   CC

      fn_del = width(1)/((W-mass(1))**2 
     &               + 0.25*width(1)*width(1))
      fn_s11 = width(2)/((W-mass(2))**2 
     &               + 0.25*width(2)*width(2))
      fn_f15 = width(3)/((W-mass(3))**2 
     &               + 0.25*width(3)*width(3))
      fn_4   = width(4)/((W-mass(4))**2 
     &               + 0.25*width(4)*width(4))

      fn_del = height(1)*fn_del
      fn_s11 = height(2)*fn_s11
      fn_f15 = height(3)*fn_f15
      fn_4   = height(4)*fn_4

c      roper_height = roper_height*(1.+q2/roper_mparm)**(-1.*roper_exp)
c      fn_roper = roper_width/((W-roper_mass)**2 
c     &               + 0.25*roper_width*roper_width) 
c      fn_roper = roper_height*fn_roper


      do i=1,4
       do j=1,4
c         fn_nr = fn_nr+
c     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))

         fn_nr = fn_nr+
     &      nr_coef(i,j)*q2**(float(j-1))*sqrt(wdif**(float(i)))
     &                  /w2temp/xb

c         write(6,*) sf

c         if(sf.EQ.2) fn_nr = fn_nr/xb
                         
       enddo
      enddo

      fn = fn_del + fn_s11 + fn_f15 + fn_4 + fn_nr
      if(sf.EQ.2) then
        fn = fn*(1.-exp(-q2/roper_mass))
      endif

c      write(6,*) "IN model:  ",sf,w2,q2,fn

      RETURN 
      END 

      subroutine christy705(W2,Q2,F1,R)
!------------------------------------------------------------------------
! subroutine to return proton structure function F1 and ratio R=sigl/sigt
!
! inputs are electron missing mass squared W2 (GeV**2) 
!            momentum trasfer squared Q2 (GeV**2)
! inputs and outputs are Real*8
! the file christy.dat is needed to use this subroutine
! Note: if W2<1.155 GeV**2, values of zero are returned (below threshold)
! Fit done by Eric Christy 7/05
! Reference this fit as ???       
!------------------------------------------------------------------------
      IMPLICIT NONE

      real*8 w2,q2,xval1(50),xvall(50),temp(4)
      real*8 mp,mp2,pi,alpha,xb,F1,FL,R
      integer i
      logical first/.true./
 
      mp = .93828
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.

      if(first) then
        first=.false.
        open(unit=15,file='f1parms.dat',status='old')
        do i=1,50
          read(15,*) temp(1)  ! par #
          read(15,*) XVAL1(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(unit=15)
        open(unit=15,file='flparms.dat',status='old')
        do i=1,50
          read(15,*) temp(1)  ! par #
          read(15,*) XVALL(i) ! starting value
          read(15,*) temp(2)   ! initial step (0 means fixed parm)
          read(15,*) temp(3) ! low limit
          read(15,*) temp(4) ! high limit 
        enddo
        close(15)
      endif

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2/(w2+q2-mp2)
      if(xb.le.0.0) return

      call resmod705(1,w2,q2,XVAL1,F1)
      call resmod705(2,w2,q2,XVALL,FL)
      if(F1.le.0.0) return

       if(F1.le.0.0) return

      F1 = F1/8.d0/pi/pi/alpha/0.389d3
      F1 = F1*(w2-mp2)
      FL = FL/8.d0/pi/pi/alpha/0.389d3
      FL = FL*(w2-mp2)*2.d0*xb
      R = FL/ (2.D0 * XB * F1)

      return
      end

CCC  Version 072205  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMOD705(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.93828
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w


      if(sf.EQ.1) then
        q2low = 0.15
      else
        q2low = 0.05
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + 0.1
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

        enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             
          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.61)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        if(sf.EQ.1) then 
         br_21_1 = 0.5
         br_21_2 = 0.5
        else
          br_21_1 = 0.985
          br_21_2 = 1.-br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

          width(6)=intwidth(6)*
     &     (xval(48)*petacm/petacmr(5)+(1.-xval(48))*ppi2cm/ppi2cmr(5))

          sig_4 =  (xval(48)*petacm+(1.-xval(48))*ppi2cm)/2./kcm/
     &           ((W2 - mass(6)**2.)**2. + (mass(6)*width(6))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

c        if(sf.EQ.2) then
c          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
c        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.(mp+mpi)**2.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 

CCC  Version 031606  -  Author:  M.E. Christy                        CCC
C*** changed to 7/7/06 version (see below)
CCC  Subroutine to get Transvese and Longitudinal eP cross sections  CCC 
CCC  from fits to L/T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


      SUBROUTINE CHRISTY31606(W2,Q2,F1,R)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval1 / 
c    &  0.12291E+01,0.15000E+01,0.15117E+01,0.16926E+01,0.16200E+01,
c    &  0.14261E+01,0.14100E+00,0.23398E+00,0.81124E-01,0.96627E-01,
c    &  0.18885E+00,0.26660E+00,0.85509E+02,0.10378E+02,0.34206E+00,
c    &  0.38546E+01,0.73111E+01,0.57629E+00,0.63691E-01,0.24072E+01,
c    &  0.19496E+02,0.38962E+01,0.25672E+00,0.29348E+01,0.10138E+02,
c    &  0.86788E+01,0.25859E+00,0.29227E+01,0.28345E+02,0.48552E+01,
c    &  0.13768E+00,0.49746E+01,0.36645E+01,0.96196E+01,0.29507E+00,
c    &  0.23934E+01,0.28438E+03,0.10536E+00,0.12390E+01,0.14593E+00,
c    &  -.24992E+03,0.67379E+00,0.23922E+01,-.24997E+00,-.18421E-02,
c    &  0.63658E-01,0.19539E+01,0.25343E+00,0.14434E+02,0.10982E+00 /
C  Simona new - h2model 7/7/06
c    & 0.12293E+01,0.15000E+01,0.15134E+01,0.16932E+01,0.16150E+01,
c    & 0.14201E+01,0.14247E+00,0.23399E+00,0.81908E-01,0.98779E-01,
c    & 0.14974E+00,0.25805E+00,0.86182E+02,0.10540E+02,0.31745E+00,
c    & 0.38879E+01,0.13214E+02,0.64008E+00,0.83161E-01,0.26896E+01,
c    & 0.18654E+02,0.45807E+01,0.24799E+00,0.30223E+01,0.10005E+02,
c    & 0.10525E+02,0.27374E+00,0.30040E+01,0.21101E+02,0.50802E+01,
c    & 0.12346E+00,0.51132E+01,0.17793E+01,0.21694E+02,0.27851E+00,
c    & 0.24741E+01,0.27502E+03,0.92225E-01,0.12263E+01,0.13867E+00,
c    & -.29351E+03,0.75890E+00,0.26699E+01,-.35054E+00,-.12799E-02,
c    & 0.73801E-01,0.20070E+01,0.51441E+00,0.30894E+02,0.94466E-01 /
C  Simona new - Christy
c    & 0.12293E+01,0.15000E+01,0.15133E+01,0.16927E+01,0.16150E+01,
c    & 0.14141E+01,0.14185E+00,0.23400E+00,0.78118E-01,0.95281E-01,
c    & 0.16498E+00,0.24073E+00,0.86154E+02,0.10245E+02,0.35896E+00,
c    & 0.38298E+01,0.12709E+02,0.67233E+00,0.85973E-01,0.25172E+01,
c    & 0.17335E+02,0.45604E+01,0.23752E+00,0.30148E+01,0.98678E+01,
c    & 0.94126E+01,0.27697E+00,0.29372E+01,0.23114E+02,0.49613E+01,
c    & 0.20701E+00,0.48640E+01,0.17441E+01,0.19868E+02,0.27348E+00,
c    & 0.24772E+01,0.27779E+03,0.89228E-01,0.12219E+01,0.14308E+00,
c    & -.29694E+03,0.75738E+00,0.26442E+01,-.35376E+00,-.11288E-02,
c    & 0.78975E-01,0.20200E+01,0.58638E+00,0.34075E+02,0.89893E-01 /
c Iteration of July 26, 2006
     & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
     & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
     & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
     & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
     & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
     & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
     & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
     & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
     & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
     & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


      data xvalL/
c    &  0.12438E+01,0.14948E+01,0.13100E+01,0.17300E+01,0.16700E+01,
c    &  0.13950E+01,0.76008E-01,0.53212E-01,0.14228E+00,0.82749E-01,
c    &  0.60000E-01,0.42000E+00,0.28298E+02,0.36911E+04,0.98237E+04,
c    &  0.11796E-03,0.63566E+01,0.26476E+05,0.36530E+05,0.33431E-03,
c    &  0.95795E+05,0.68291E+04,0.35634E+04,0.90878E+05,0.12332E+03,
c    &  0.73314E+06,0.26074E+06,0.74001E+02,0.00000E+00,0.92972E+04,
c    &  0.87089E+04,0.46171E+05,0.10478E+02,0.73525E+04,0.95318E+04,
c    &  0.29686E+00,0.69642E+03,0.13502E+01,0.00000E+00,0.34422E+01,
c    &  0.19016E+01,0.20149E+00,0.62176E+03,0.61984E+01,0.16046E+00,
c    &  0.16765E+01,0.00000E+00,0.40092E-11,0.71919E+00,0.33674E+00 /
C  Simona new - h2model 7/7/06
c    & 0.12437E+01,0.14943E+01,0.12808E+01,0.17330E+01,0.16700E+01,
c    & 0.13921E+01,0.71753E-01,0.45000E-01,0.57526E-01,0.73573E-01,
c    & 0.60000E-01,0.47000E+00,0.28690E+02,0.73923E+04,0.20383E+05,
c    & 0.00000E+00,0.53160E+01,0.49592E+05,0.64624E+05,0.00000E+00,
c    & 0.12151E+05,0.10000E+05,0.31121E+04,0.81442E+05,0.41594E+04,
c    & 0.38606E+06,0.96177E+05,0.42609E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17969E+02,0.69923E+04,0.75455E+04,
c    & 0.79081E+00,0.69950E+03,0.13591E+01,0.00000E+00,0.34511E+01,
c    & 0.19061E+01,0.22893E+00,0.10577E+04,0.79123E+01,0.12661E+00,
c    & 0.15646E+01,0.00000E+00,0.00000E+00,0.66613E+00,0.38768E+00 /
C  Simona new - Christy
c    & 0.12441E+01,0.14943E+01,0.12820E+01,0.17304E+01,0.16700E+01,
c    & 0.13955E+01,0.73364E-01,0.45064E-01,0.68170E-01,0.76574E-01,
c    & 0.60000E-01,0.46067E+00,0.29174E+02,0.73983E+04,0.20326E+05,
c    & 0.00000E+00,0.53291E+01,0.49445E+05,0.64863E+05,0.00000E+00,
c    & 0.13790E+05,0.10000E+05,0.37341E+04,0.61390E+05,0.44336E+04,
c    & 0.37016E+06,0.10266E+06,0.39083E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.16753E+02,0.70252E+04,0.75339E+04,
c    & 0.77189E+00,0.70838E+03,0.13292E+01,0.00000E+00,0.34217E+01,
c    & 0.19038E+01,0.21016E+00,0.98783E+03,0.80689E+01,0.12300E+00,
c    & 0.15563E+01,0.00000E+00,0.00000E+00,0.65655E+00,0.36031E+00 /
C  Simona - July 28, 2006
     & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
     & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
     & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
     & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
     & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
     & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
     & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
     & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
     & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
     & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmod316(1,w2,q2,xval1,sigT)
      call resmod316(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

   
      end

         SUBROUTINE H2MODEL(QQ4,WW4,W1,W2)

********************************************************************************
*
* This subroutine calculates model cross sections for inclusive electron-proton
* scattering in the resonance region. The cross section is returned in 
* nanobarns/sr GeV. This fit is a modified version of Linda Stuart's 8/91 fit. 
* One major difference is that the coefficients represent fit results from a 
* substantially expanded data set, including all inclusive SLAC data in the 
* range 1.1 < W^2 < 5. There are other differences; for a complete discussion, 
* see Keppel's Ph.D. thesis. 2/94 CEK
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
* SIG1     = Cross section in nb/sr/GeV**2 (DSIG/DOMEGA/DW**2)
* SIG2     = Cross section in nb/sr/GeV (DSIG/DOMEGA/DEP)
* SIG_NRES = Non-resonant contribution to SIG1.
* SIG_RES  = Resonant contribution to SIG1.
* SIG_NR   = Non-resonant contribution to SIG2.
* SIG_R    = Resonant contribution to SIG2.
* SIGroper = Possible Roper contribution to SIG2.
* goroper  = Logical variable, set true if including possible Roper strength
*
* SIG1, SIG2, SIG_NRES, and SIG_RES are 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
        IMPLICIT NONE

        logical goroper 
        logical goodfit 
        INTEGER I
        REAL*4  QQ4,WW4,W1,W2
        REAL*8  SIN2, SIG1(2), SIG2(2), SIG_RES(2), SIG_R(2), 
     >          SIG_RES1(2), SIG_RES2(2), SIG_NRES(2), SIG_NR(2), 
     >          DIPOLE, COS2, TAU, EPS, K, DW2DEP, DEPCONV, 
     >          PI, AM, ALPHA,NU,CONV,SIGT,SIGL,qq,ww,
     >          RADCON, R_NRES, sigroper(2), x, r, dr,fac,rlog
       
        qq = qq4
        ww = ww4

        RADCON = 3.141592654/180.
        PI = 3.14159265
        AM = .9382727
        ALPHA = 0.00729735
        CONV = 0.0025767

C        SIN2 = SIN(TH*RADCON/2.0)*SIN(TH*RADCON/2.0)
C        COS2 = 1 - SIN2
C        Q2 = 4.0*E*EP*SIN2
C        W2 = AM*AM + 2.0*AM*(E - EP) - 4.0*E*EP*SIN2

         goroper = .true.

C        write(6,*) q2,w2

        IF(WW.LT.1.15) THEN            ! Below pion threshold
          DO I = 1,2
            SIG1(I) = 0.0
            SIG2(I) = 0.0
            SIG_NRES(I) = 0.0
            SIG_RES(I) = 0.0
           SIG_NR(I) = 0.0
            SIG_R(I) = 0.0
            sigroper(I) = 0.0
          ENDDO
          RETURN
        ENDIF

        K = (WW - AM*AM)/(2.0*AM)
        NU = K + QQ/(2.0*AM)
        TAU = NU*NU/QQ
        x = qq/2./AM/nu
C        EPS = 1.0/(1.0 + 2.0*(1.0 + TAU)*SIN2/COS2)
        DIPOLE = 1.0/(1.0 + QQ/0.71)**2
c        DEPCONV = ALPHA*K*EP/(4.0*PI*PI*Q2*E)*(2.0/(1.0 - EPS))*1000.
c        DW2DEP = 2.0*AM + 4.0*E*SIN2

! H2MOD_FIT returns cross sections in units of microbarns/(dipole FF)**2
        CALL H2MODEL_FIT(QQ,WW,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)
c        call r1990(x,qq,r,dr,goodfit)
c        if (goodfit) then
c         R_NRES = r
c        else 
         R_NRES = 0.25/SQRT(QQ)
c        endif

        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
C        SIG_NRES(2) = SIG_NRES(2)*DIPOLE*DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1))
     >                 *DIPOLE*DIPOLE
C        SIG_RES(2)  = SQRT(SIG_RES1(2)**2 + SIG_RES2(2)**2 + 
C     >                 sigroper(2)**2)*
C     >                 DIPOLE*DIPOLE
C        SIG1(1) = SIG_NRES(1)*(1.0 + EPS*R_NRES) + SIG_RES(1)
C        SIG1(2) = SQRT( (SIG_NRES(2)*(1.0 + EPS*R_NRES))**2 +
C     >                   SIG_RES(2)**2)
C        SIG2(1) = SIG1(1)*DW2DEP
C        SIG2(2) = SIG1(2)*DW2DEP
C        sig_nr(1) = sig_nres(1)*dw2dep
C        sig_nr(2) = sig_nres(2)*dw2dep
C        sig_r(1) = sig_res(1)*dw2dep
C        sig_r(2) = sig_res(2)*dw2dep
        
C        sige = sig2(1)

         SIGT = SIG_NRES(1)+SIG_RES(1)
c         sigl = r*sigt  
        SIGL = R_NRES*SIG_NRES(1)
         W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
         W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END

      SUBROUTINE H2MODEL_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J, KK
      REAL*8 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     WR, KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, GAM, 
     >     GR, RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER, ERRCHECK

      REAL*8 PI, ALPHA, AM, MPPI, MDELTA, MPI, GAM1WID, GAM2WID, MASS1,
     >     ROPERWID, MASSROPER, MASS2, DELWID, FIT(30), SQRTWDIF,
     >     XR, MQDEP

       REAL*8 RCOEF(27), COEF(22), RERR1(200),rerr2(125), ERR1(200), 
     >        ERR2(53)
    
      LOGICAL FIRST

      FIRST = .TRUE.

c      pi = 3.14159265
c      alpha = 0.00729735
c      am = .9382727
c      MPPI = 1.07783
c      MDELTA = 1.2340 
c      MPI = 0.13957 
c      GAM1WID = 0.0800 
c      GAM2WID = 0.0900
c      MASS1 = 1.5045 
c      ROPERWID = 0.0500 
c      MASSROPER = 1.4000 
c      MASS2 = 1.6850
c      DELWID = 0.1200
c      XR = 0.1800 
c      MQDEP = 3.40 
      
      pi = 3.14159265
      alpha = 0.00729735
      am = .9382727
      MPPI = 1.07783
      MDELTA = 1.229 
      MPI = 0.13957 
      GAM1WID = 0.080
      GAM2WID = 0.090
      MASS1 = 1.5062 
      ROPERWID = 0.0500 
      MASSROPER = 1.4000 
      MASS2 = 1.6810
      DELWID = 0.120
      XR = 0.1800 
      MQDEP = 3.40 

      DATA RCOEF/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >   4.6273E-01,   1.7844E-01,   -1.5416E+02, -1.4891E+02,
     >   2.4102E+02,   2.5823E+00,    7.1004E+00, -8.9771E+00,
     >   1.3744E+00,  -1.2085E+00,    1.1218E-01/

      DATA COEF/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/
 
c      DATA RERR1/
c     >  2.6120E+02,-9.4211E+02, 4.0844E+03, 7.4994E+02,-3.5317E+03,
c     >  3.1703E+03, 2.1238E-01,-6.1568E-01, 4.1479E-01, 1.9720E-02,
c     >  1.2891E-01,-4.1615E+00, 4.7246E+00, 2.8090E-03, 6.1657E-02,
c     >  1.3120E+00,-9.4379E+00, 9.0902E+00,-1.3636E-03, 2.8054E-02,
c     >  9.6123E-02,-1.1465E+03, 3.9099E+03,-3.0097E+03,-1.0604E+01,
c     > -1.2214E+00,-8.3549E-01, 1.3696E+04, 3.9085E+03,-1.5369E+04,
c     >  1.2745E+04, 2.9942E+01, 7.7268E+00, 1.0113E+01,-4.3868E+04,
c     >  1.5709E+05,-3.0207E+03, 1.2809E+04,-1.1075E+04,-2.0442E+01,
c     > -7.5843E+00,-1.0773E+01, 3.2597E+04,-1.2457E+05, 1.0283E+05,
c     > -1.6960E-01, 5.9410E-01,-4.5486E-01,-1.0715E-02,-2.6512E-03,
c     >  1.0153E-03, 6.4074E+00,-1.9189E+01, 1.3612E+01, 9.2486E-03,
c     >  2.7904E-01, 6.3576E+00,-7.8552E+00, 1.5302E-02,-1.1506E-01,
c     > -4.7552E-02,-1.0171E+01,-1.5884E+00, 1.6223E+01,-1.1379E-04,
c     >  4.9212E-01,-2.4354E+00, 1.7921E+01,-1.7223E+01, 4.0778E-03,
c     > -4.5558E-02,-1.8539E-01, 7.9930E+00,-7.1588E+01, 7.1512E+01,
c     > -2.1529E-03, 1.8337E-01, 7.7590E-01, 7.3007E+02,-2.5219E+03,
c     >  1.9547E+03, 6.1102E+00, 1.2970E+00,-1.3084E+00,-9.4932E+03,
c     >  3.1542E+04,-2.3894E+04,-5.9583E+00, 8.1005E-02, 3.6885E-01,
c     >  9.5708E+03,-2.4911E+03, 9.4342E+03,-7.7120E+03,-1.8608E+01,
c     > -1.1065E+00, 6.5015E+00, 3.1755E+04,-1.1529E+05, 9.1964E+04,
c     >  1.8347E+01,-2.5899E+00, 7.1169E-01,-3.2268E+04, 1.1891E+05,
c     >  1.9339E+03,-7.7737E+03, 6.6128E+03, 1.3392E+01,-7.3587E-02,
c     > -4.9353E+00,-2.4223E+04, 9.2768E+04,-7.6712E+04,-1.3210E+01,
c     >  1.2513E+00,-4.5156E+00, 2.4541E+04,-9.5131E+04, 7.8848E+04,
c     >  1.9952E-02,-7.1332E-02, 5.5522E-02, 9.8804E-04, 2.3682E-04,
c     > -7.9762E-05,-6.3638E-01, 1.9492E+00,-1.4036E+00,-9.9312E-04,
c     > -7.8634E-05, 8.2617E-05, 6.8002E-01,-2.1138E+00, 1.5308E+00,
c     >  1.3008E-04,-1.0389E+02, 3.5942E+02,-2.7883E+02,-6.0671E-01,
c     > -1.3016E-01, 1.4621E-01, 1.2841E+03,-4.3361E+03, 3.3132E+03,
c     >  7.0701E-01, 1.2805E-01, 1.3355E-01,-1.4645E+03, 4.9522E+03,
c     > -3.7686E+03,-1.0047E-01, 2.7406E+02, 3.5483E+02,-1.3433E+03,
c     >  1.0978E+03, 1.9033E+00, 5.3726E-02,-8.1621E-01,-4.3612E+03,
c     >  1.6110E+04,-1.2957E+04,-2.2247E+00,-2.1299E-01,-5.8178E-01,
c     >  4.9755E+03,-1.8393E+04, 1.4724E+04, 3.1774E-01,-9.2555E+02,
c     >  3.4086E+03,-2.7508E+02, 1.1025E+03,-9.3653E+02,-1.4100E+00,
c     >  7.3163E-02, 6.6492E-01, 3.3590E+03,-1.3073E+04, 1.0893E+04,
c     >  1.6311E+00, 2.4826E-01, 8.3308E-01,-3.7999E+03, 1.4772E+04,
c     > -1.2252E+04,-2.3255E-01, 7.0167E+02,-2.7162E+03, 2.2434E+03,
c     >  3.0688E+00,-1.0328E+01, 7.8828E+00, 3.6601E-03, 1.3367E-03,
c     > -2.9672E-03,-3.2441E+01, 1.0979E+02,-8.3795E+01,-6.6345E-03/

c      DATA  rerr2/
c     > 3.7074E-02,
c     >-5.7300E-02, 1.5212E-02, 4.5952E-04,
c     > 1.1568E-04,-2.9315E-04,-4.6018E-01,
c     > 9.3624E-01,-4.5908E-01,-6.2914E-05,
c     > 1.1699E-03, 2.0141E-03, 6.9968E-02,
c     >-1.9348E-01, 1.2176E-01, 5.4214E-07,
c     > 1.3845E-04, 2.5311E-03,-2.5396E-03,
c     >-1.2757E-04, 2.4971E-04,-1.2737E-04,
c     > 7.2023E-03,-4.1404E-03, 4.6704E-04,
c     > -4.6388E-03,-5.2545E-03, 4.0159E+01,-1.3481E+02, 1.0186E+02,
c     >  1.1796E-03,-9.1088E+00, 3.0200E+01,-2.2552E+01, 4.3562E-01,
c     > -1.0404E+01, 3.8414E+01,-3.0978E+01,-1.4730E-02, 4.6327E-03,
c     >  1.9716E-02, 1.1236E+02,-4.1952E+02, 3.3862E+02, 2.4150E-02,
c     >  1.1098E-02, 2.0122E-02,-1.3812E+02, 5.1058E+02,-4.0773E+02,
c     > -4.1791E-03, 3.0702E+01,-1.1132E+02, 8.7622E+01,-1.4199E+00,
c     >  5.0230E+00, 8.0171E+00,-3.1384E+01, 2.6350E+01, 1.3147E-02,
c     > -6.1508E-03,-1.6808E-02,-8.7538E+01, 3.4530E+02,-2.8922E+02,
c     > -1.9581E-02,-1.0895E-02,-2.4705E-02, 1.0611E+02,-4.1369E+02,
c     >  3.4296E+02, 3.2847E-03,-2.3191E+01, 8.8502E+01,-7.2288E+01,
c     >  1.0469E+00,-3.8873E+00, 3.1142E+00,
c     > 1.1348E+00,-1.7657E+00, 4.7686E-01,
c     > 1.6653E-02, 4.3488E-04,-7.5168E-03,
c     >-1.6696E+01, 3.4692E+01,-1.7470E+01,
c     >-4.9697E-03, 4.4232E-02, 5.7617E-02,
c     > 5.7800E+00,-1.3886E+01, 7.9819E+00,
c     > 3.4744E-04,-5.4411E-01, 1.2683E+00,
c     >-7.0771E-01, 1.1282E-02,-2.4800E-02,
c     > 1.2909E-02, 1.5171E-01,-6.0417E-01,
c     > 7.7405E-01,-5.8981E-02,-5.8502E-03,
c     > 8.8611E-04, 5.8326E-03, 6.5418E+00,
c     >-1.2978E+01, 6.1069E+00, 1.2462E-03,
c     >-1.8442E-02,-2.7954E-02,-1.8335E+00,
c     > 4.3674E+00,-2.4393E+00,-6.2354E-05,
c     > 1.4746E-01,-3.4127E-01, 1.8285E-01,
c     >-3.0479E-03, 6.8138E-03,-3.4673E-03,
c     >-7.5270E-02, 4.0914E-02/


c      DATA ERR1/
c     >  3.7797E+02,-1.2732E+03, 4.8470E+03, 9.7589E+02,-3.9592E+03,
c     >  3.3447E+03, 1.9629E-01,-4.2402E-01, 1.9757E-01, 3.0613E-02,
c     > -4.0257E-01,-2.0922E+00, 3.0126E+00, 3.8385E-03, 7.3553E-02,
c     >  1.4084E+00,-8.4718E+00, 7.8586E+00,-1.6484E-03, 2.2185E-02,
c     >  7.4896E-02,-1.5627E+03, 5.0106E+03,-3.7125E+03,-1.1701E+01,
c     > -6.9186E-01,-1.4263E+00, 1.5792E+04, 5.0288E+03,-1.7793E+04,
c     >  1.3974E+04, 3.1643E+01, 5.0040E+00, 9.9958E+00,-4.8540E+04,
c     >  1.6247E+05,-3.7498E+03, 1.4066E+04,-1.1461E+04,-2.0806E+01,
c     > -5.0428E+00,-9.7813E+00, 3.5056E+04,-1.2382E+05, 9.7850E+04,
c     > -2.0038E-01, 5.9769E-01,-4.0397E-01,-1.5776E-02,-3.7509E-03,
c     >  5.7496E-04, 7.2218E+00,-2.0335E+01, 1.3722E+01, 1.2562E-02,
c     >  1.4708E+00, 1.8510E+00,-4.1856E+00, 1.9572E-02,-1.3469E-01,
c     > -3.7791E-02,-1.5215E+01, 1.8843E+01,-9.9384E-01, 5.4133E-04,
c     >  5.6775E-01,-2.4158E+00, 1.5245E+01,-1.4180E+01, 5.3668E-03,
c     > -3.5419E-02,-1.4360E-01, 7.8707E+00,-5.7677E+01, 5.5406E+01,
c     > -7.5727E-04, 1.4127E-01, 5.8964E-01, 1.0277E+03,-3.3407E+03,
c     >  2.4943E+03, 6.1372E+00, 2.0731E+00,-1.0628E-01,-1.1445E+04,
c     >  3.6033E+04,-2.6376E+04,-6.4849E+00,-1.5437E+00,-3.1093E+00,
c     >  1.1966E+04,-3.3062E+03, 1.1473E+04,-8.9323E+03,-1.7658E+01,
c     > -3.0298E+00, 2.4862E+00, 3.6140E+04,-1.2237E+05, 9.3797E+04,
c     >  1.8377E+01, 2.4649E-01, 9.5713E+00,-3.7362E+04, 1.2613E+05,
c     >  2.4733E+03,-8.9836E+03, 7.2301E+03, 1.2133E+01, 1.0120E+00,
c     > -2.0972E+00,-2.6581E+04, 9.4364E+04,-7.4804E+04,-1.2397E+01,
c     >  5.8276E-01,-9.1893E+00, 2.7145E+04,-9.6250E+04, 7.6086E+04,
c     >  2.4070E-02,-7.3772E-02, 5.1165E-02, 1.4597E-03, 3.3977E-04,
c     > -2.6275E-05,-7.2542E-01, 2.0676E+00,-1.4052E+00,-1.3577E-03,
c     > -1.4477E-04,-8.5451E-05, 7.4811E-01,-2.1217E+00, 1.4288E+00,
c     >  1.7439E-04,-1.6022E+02, 5.2231E+02,-3.9172E+02,-4.1771E-01,
c     > -2.3133E-01,-1.9119E-02, 1.6931E+03,-5.4146E+03, 4.0099E+03,
c     >  6.5228E-01, 4.5766E-01, 6.7254E-01,-2.0266E+03, 6.3551E+03,
c     > -4.6404E+03,-9.4689E-02, 4.2768E+02, 5.1531E+02,-1.7829E+03,
c     >  1.3890E+03, 1.1798E+00, 3.1335E-01,-2.5902E-01,-5.3955E+03,
c     >  1.8502E+04,-1.4311E+04,-1.8045E+00,-9.6753E-01,-2.0260E+00,
c     >  6.3626E+03,-2.1445E+04, 1.6387E+04, 2.6350E-01,-1.3456E+03,
c     >  4.5055E+03,-3.8598E+02, 1.3911E+03,-1.1170E+03,-7.9328E-01,
c     > -7.6840E-02, 2.5967E-01, 4.0005E+03,-1.4347E+04, 1.1455E+04,
c     >  1.1795E+00, 6.2629E-01, 1.6961E+00,-4.6485E+03, 1.6399E+04,
c     > -1.2954E+04,-1.7187E-01, 9.8638E+02,-3.4363E+03, 2.7002E+03,
c     >  6.0266E+00,-1.9528E+01, 1.4686E+01,-1.7956E-02, 3.3364E-03,
c     >  1.2080E-03,-5.5018E+01, 1.7933E+02,-1.3517E+02, 7.9955E-03/


c       DATA ERR2/
c     > -2.1546E-02,-2.3493E-02, 7.4315E+01,-2.3518E+02, 1.7398E+02,
c     > -6.4429E-04,-1.9950E+01, 6.3147E+01,-4.6881E+01, 1.2816E+00,
c     > -1.9366E+01, 6.5755E+01,-5.0971E+01, 5.7005E-02, 3.3439E-04,
c     >  5.5786E-03, 1.7715E+02,-6.1369E+02, 4.7999E+02,-2.9558E-02,
c     >  5.5461E-02, 7.1075E-02,-2.3560E+02, 7.9154E+02,-6.0792E+02,
c     >  2.7242E-03, 6.3265E+01,-2.0981E+02, 1.6050E+02,-4.0749E+00,
c     >  1.3388E+01, 1.4562E+01,-5.1058E+01, 4.0557E+01,-4.3474E-02,
c     > -4.4868E-03,-6.3041E-03,-1.3274E+02, 4.7814E+02,-3.8441E+02,
c     >  2.5677E-02,-3.8538E-02,-5.8204E-02, 1.7424E+02,-6.0799E+02,
c     >  4.8014E+02,-2.6425E-03,-4.6992E+01, 1.6058E+02,-1.2570E+02,
c     >  3.0554E+00,-1.0258E+01, 7.9929E+00/


! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
!Define error matrix:
c        KK = 0
c        if (goroper) then 
c          DO J = 1,25
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) RERRMAT(I,J) = RERR1(KK)
c              if (kK.le.325.and.kK.gt.200) RERRMAT(I,J)=RERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c       if (.not.goroper) then  
c          DO J = 1,22
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) ERRMAT(I,J) = ERR1(KK)
c              if (kK.le.253.and.kK.gt.200) ERRMAT(I,J)=ERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c        if (goroper) then
c          DO J = 1,25
c              DO I = J+1,25
c                RERRMAT(I,J) = RERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
c        if (.not.goroper) then
c          DO J = 1,22
c              DO I = J+1,22
c                ERRMAT(I,J) = ERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(25) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
        FIT(I + 12) = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 18)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 21)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO
      if (goroper) FIT(26)  = FIT(25)/sqrt(Q2)
      if (goroper) FIT(27)  = FIT(25)/Q2

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,27
c          RSIGERR(J) = 0.0
c          DO I = 1,25
c            RSIGERR(J) = RSIGERR(J) + FIT(J)*FIT(I)*RERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*RERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12.OR.J.EQ.17
     > .OR.J.EQ.18 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + RSIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
c             SIG_RES1(2) = SIG_RES1(2) + RSIGERR(J)
          elseIF(j.ge.25.and.j.le.27) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
c            SIGroper(2) = SIGroper(2) + RSIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + RSIGERR(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
c          SIGERR(J) = 0.0
c          DO I = 1,22
c            SIGERR(J) = SIGERR(J) + FIT(J)*FIT(I)*ERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*ERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + SIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
c            SIG_RES1(2) = SIG_RES1(2) + SIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + SIGERR(J)
          ENDIF
        ENDDO
      endif

! ERRCHECK should agree with SIGTOTER.
C      ERRCHECK = SQRT(ABS(SIG_RES2(2) + SIG_RES1(2) + SIG_NRES(2)))
C      SIGTOTER = SQRT(SIGTOTER)
c      SIG_RES2(2) = SQRT(ABS(SIG_RES2(2)))
c      SIG_RES1(2) = SQRT(ABS(SIG_RES1(2)))
c      SIG_NRES(2) = SQRT(ABS(SIG_NRES(2)))

      RETURN
      END
CCC  Version 050521  -  Author:  M.E. Christy   modified by S. Malace CCC
CCC  Subroutine to get Transvese eD cross sections  CCC 
CCC  from fits to T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


C S. Malace 7/22/06 to get D2 

      SUBROUTINE rescssim(W2,Q2,F1)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval11(50),xval12(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1
      integer i,npts,sf
 
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11872E+01,0.15218E+01,0.14125E+01,0.21081E+01,0.17094E+01,
c    & 0.28094E+01,0.27909E+00,0.13229E+00,0.43753E+00,0.54004E+00,
c    & 0.30149E+00,0.97815E+00,0.491276E+04,0.188911E+04,0.352732E+04,
c    & -.29889E+03,0.253598E+04,0.96323E+03,0.881774E+05,0.344864E+05,
c    & 0.28063E+02,0.633321E+05,0.761637E+05,-.473462E+04,0.90734E+02,
c    & 0.8165041E+06,0.264528E+06,0.808716E+05,0.458931E+05,0.88573E+03,
c    & 0.12218E+02,0.37465E+01,0.90133E+03,0.1007504E+05,0.4736814E+05,
c    & -.463262E+04,0.356829E+03,0.17208E+01,0.13029E+01,0.65986E-01,
c    & 0.6864799E+05,0.375928E+02,0.342095E+01,-.91793E+00,0.23611E-02,
c    & 0.27051E+00,0.136559E+01,0.65689E+00,0.171501E+03,0.26701E+01 /
c
c
c     data xval12/
c
c
C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11996E+01,0.153641E+01,0.142073E+01,0.212016E+01,0.171082E+01,
c    & 0.287999E+01,0.29076E+00,0.192957E+00,0.439143E+00,0.557558E+00,
c    & 0.326104E+00,0.97999E+00,0.45174E+04,0.191460E+04,0.353009E+04,
c    & -.31843E+03,0.600311E+04,0.87710E+03,0.8694928E+05,0.3897957E+05,
c    & 0.28447E+02,0.6324006E+05,0.7590074E+05,-.48051E+04,0.99998E+02,
c    & 0.8355478E+06,0.2756187E+06,0.845187E+05,0.490728E+05,0.9423E+03,
c    & 0.13145E+02,0.36846E+01,0.105307E+04,0.96985E+04,0.472598E+05,
c    & -.46531E+04,0.35792E+03,0.1734E+01,0.13036E+01,0.65573E-01,
c    & 0.908018E+05,0.395258E+02,0.341936E+01,-.90659E+00,0.24026E-02,
c    & 0.27467E+00,0.136335E+01,0.66170E+00,0.16993E+03,0.26701E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

cc Iteration2 follows: out in simonaD2SimFig3.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.119098E+01,0.153249E+01,0.142649E+01,0.213592E+01,0.17119E+01,
c    & 0.297601E+01,0.271804E+00,0.173939E+00,0.394625E+00,0.568884E+00,
c    & 0.296322E+00,0.979981E+00,0.49315E+04,0.26699E+04,0.199866E+04,
c    & -.29483E+03,0.99585E+04,0.177447E+04,0.995236E+05,0.376651E+05,
c    & 0.32501E+02,0.631396E+05,0.772169E+05,-.457579E+04,0.989009E+02,
c    & 0.8603541E+06,0.2483892E+06,0.8244512E+05,0.47721E+05,0.9463E+03,
c    & 0.14093E+02,0.39168E+01,0.12428E+04,0.1017448E+05,0.476793E+05,
c    & -.46163E+04,0.358308E+03,0.17098E+01,0.13262E+01,0.65102E-01,
c    & 0.942561E+05,0.65499E+02,0.30852E+01,-.954945E+00,0.20332E-02,
c    & 0.29676E+00,0.135420E+01,0.68348E+00,0.16509E+03,0.2669E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c Following is for simonaD2SimFit4.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11888E+01,0.15245E+01,0.142628E+01,0.213139E+01,0.17127E+01,
c    & 0.29098E+01,0.30225E+00,0.141739E+00,0.41125E+00,0.550194E+00,
c    & 0.29422E+00,0.97980E+00,0.456537E+04,0.181743E+04,0.25647E+04,
c    & -.36697E+03,0.24858E+04,0.78326E+03,0.952766E+05,0.346543E+05,
c    & 0.30294E+02,0.664647E+05,0.787561E+05,-.43423E+04,0.9999E+02,
c    & 0.9195253E+06,0.3002124E+06,0.871368E+05,0.51180E+05,0.9766E+03,
c    & 0.1320E+02,0.3838E+01,0.10652E+04,0.100941E+05,0.4728109E+05,
c    & -.460501E+04,0.36059E+03,0.17467E+01,0.130749E+01,0.67122E-01,
c    & 0.965392E+05,0.36514E+02,0.336334E+01,-.86248E+00,0.21560E-02,
c    & 0.27810E+00,0.135956E+01,0.657735E+00,0.16344E+03,0.26694E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c following is iteration 3 starting with F2allm. /rc/simonaD2SimFit5.out
c     data xval11 / 

C  Simona D2 fit to (simona data) + (ioana-delta)
c    & 0.11871E+01,0.153116E+01,0.143781E+01,0.21428E+01,0.171332E+01,
c    & 0.29792E+01,0.26959E+00,0.171029E+00,0.39371E+00,0.570214E+00,
c    & 0.29644E+00,0.973246E+00,0.503525E+04,0.263216E+04,0.19497E+04,
c    & -.30213E+03,0.950003E+04,0.17909E+04,0.9990789E+05,0.374659E+05,
c    & 0.32625E+02,0.633281E+05,0.775036E+05,-.450712E+04,0.99982E+02,
c    &0.8657796E+06,0.2497579E+06,0.82912E+05,0.481666E+05,0.95763E+03,
c    & 0.14368E+02,0.39312E+01,0.13487E+04,0.101466E+05,0.476961E+05,
c    & -.46158E+04,0.35876E+03,0.17086E+01,0.132586E+01,0.64952E-01,
c    & 0.927868E+05,0.65323E+02,0.30805E+01,-.956083E+00,0.211227E-02,
c    & 0.29716E+00,0.135334E+01,0.68561E+00,0.16425E+03,0.26698E+01 /


c     data xval12/


C  Simona D2 fit to ioana data
c    & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
c    & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
c    & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
c    & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
c    & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
c    & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
c    & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
c    & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
c    & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
c    & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

c following is iteration 3 starting with Bodek. /rc/simonaD2SimFit6.out
C  Simona D2 fit to (simona data) + (ioana-delta)
      data xval11/
     & 0.118651E+01,0.15247E+01,0.14354E+01,0.213435E+01,0.17131E+01,
     & 0.29094E+01,0.30152E+00,0.14436E+00,0.400969E+00,0.55492E+00,
     & 0.29446E+00,0.97992E+00,0.460104E+04,0.17992E+04,0.25114E+04,
     & -.37399E+03,0.24733E+04,0.75449E+03,0.8939804E+05,0.3523602E+05,
     & 0.30091E+02,0.665424E+05,0.790445E+05,-.42655E+04,0.99998E+02,
     & 0.9217486E+06,0.2988304E+06,0.869219E+05,0.512361E+05,0.9782E+03,
     & 0.13294E+02,0.38278E+01,0.10849E+04,0.100891E+05,0.4726618E+05,
     & -.46058E+04,0.36094E+03,0.17463E+01,0.130750E+01,0.670575E-01,
     & 0.969567E+05,0.361745E+02,0.335147E+01,-.86421E+00,0.21919E-02,
     & 0.28124E+00,0.135992E+01,0.65969E+00,0.16259E+03,0.26693E+01 /


      data xval12/


C  Simona D2 fit to ioana data
     & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
     & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
     & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
     & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
     & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
     & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
     & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
     & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
     & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
     & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

       if(q2.lt.4.5.and.w2.lt.2.and.q2.gt.1.7) then

         call resmodsim(1,w2,q2,xval12,sigt)

       else

         call resmodsim(1,w2,q2,xval11,sigt)

      endif  

      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3

      return
      end

CCC  Version 061105  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

C THIS IS FOR D2 from S. MALACE 7/22/07
             
      SUBROUTINE RESMODSIM(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.9382727
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (0.937272 + 0.136)
      wr = wdif/w


c      if(q2.LT.0.3.AND.sf.EQ.1) then
      if(q2.LT.0.15) then
        lowq2 = .true.
        lmax = 2
      endif

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = 0.15
        elseif(l.EQ.2.AND.lowq2) then
          q2 = 0.25
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

       enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.71)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        br_21_1 = 0.5
        br_21_2 = 0.5
        if(sf.EQ.2) then
          br_21_1 = xval(48)
          br_21_2 = 1.- br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +(nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.+nr_coef(i,3)*wdif**(float(3*i-1)/2)
     &       *(1.-xth(1))**3.)/(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

        if(sf.EQ.2) then
c          sig_res = sig_res + sig_4L
          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.1.16.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(0.15-q2)
        else
          slope = sig/0.15
          sig = sig - slope*(0.15-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 

      SUBROUTINE CHRISTYX(W2,Q2,F1,R,sigt,sigl)
      
      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

c     data xval1 / 
c Iteration of July 26, 2006
c    & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
c    & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
c    & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
c    & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
c    & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
c    & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
c    & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
c    & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
c    & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
c    & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


c     data xvalL/
C  Simona - July 28, 2006
c    & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
c    & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
c    & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
c    & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
c    & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
c    & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
c    & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
c    & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
c    & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

c copied from 806 version
c     data xval1 / 
c    & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
c    & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
c    & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
c    & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
c    & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
c    & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
c    & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
c    & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
c    & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
c    & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


c     data xvalL/
c    & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
c    & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
c    & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
c    & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
c    & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
c    & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
c    & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
c    & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
c    & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
c    & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

c from original version
      data xval1 / 
     &  0.12291E+01,0.15000E+01,0.15117E+01,0.16926E+01,0.16200E+01,
     &  0.14261E+01,0.14100E+00,0.23398E+00,0.81124E-01,0.96627E-01,
     &  0.18885E+00,0.26660E+00,0.85509E+02,0.10378E+02,0.34206E+00,
     &  0.38546E+01,0.73111E+01,0.57629E+00,0.63691E-01,0.24072E+01,
     &  0.19496E+02,0.38962E+01,0.25672E+00,0.29348E+01,0.10138E+02,
     &  0.86788E+01,0.25859E+00,0.29227E+01,0.28345E+02,0.48552E+01,
     &  0.13768E+00,0.49746E+01,0.36645E+01,0.96196E+01,0.29507E+00,
     &  0.23934E+01,0.28438E+03,0.10536E+00,0.12390E+01,0.14593E+00,
     &  -.24992E+03,0.67379E+00,0.23922E+01,-.24997E+00,-.18421E-02,
     &  0.63658E-01,0.19539E+01,0.25343E+00,0.14434E+02,0.10982E+00 /


      data xvalL/
     &  0.12438E+01,0.14948E+01,0.13100E+01,0.17300E+01,0.16700E+01,
     &  0.13950E+01,0.76008E-01,0.53212E-01,0.14228E+00,0.82749E-01,
     &  0.60000E-01,0.42000E+00,0.28298E+02,0.36911E+04,0.98237E+04,
     &  0.11796E-03,0.63566E+01,0.26476E+05,0.36530E+05,0.33431E-03,
     &  0.95795E+05,0.68291E+04,0.35634E+04,0.90878E+05,0.12332E+03,
     &  0.73314E+06,0.26074E+06,0.74001E+02,0.00000E+00,0.92972E+04,
     &  0.87089E+04,0.46171E+05,0.10478E+02,0.73525E+04,0.95318E+04,
     &  0.29686E+00,0.69642E+03,0.13502E+01,0.00000E+00,0.34422E+01,
     &  0.19016E+01,0.20149E+00,0.62176E+03,0.61984E+01,0.16046E+00,
     &  0.16765E+01,0.00000E+00,0.40092E-11,0.71919E+00,0.33674E+00 /

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmodX(1,w2,q2,xval1,sigT)
      call resmodX(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

   
      end

CCC  Version 031606  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMODX(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low,dq2,t,xpr
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2
      common/tst1/sigr,sig_nr

      lowq2 = .false.
      lmax = 1
      q2temp = q2
      dq2 = 0.05

      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w

      br(1,1) = 1.0     !!! single pion branching ratios
      br(2,1) = 0.5
      br(3,1) = 0.65
      br(4,1) = 0.65
      br(5,1) = 0.4
      br(6,1) = 0.65
      br(7,1) = 0.6

      if(sf.EQ.2) then 
        br(6,1) = xval(48)
        br(2,1) = xval(49)
      endif 

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  ? 4th resonance region

      do i=1,7
        x0(i) = 0.165
      enddo
      x0(4) = 0.6

      do i=1,7
        br(i,2) = 1.-br(i,1)
      enddo
    

      if(sf.EQ.1) then
        q2low = 0.00
      else
        q2low = 0.1
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + dq2
        endif

        dip = 1./(1.+q2/0.71)**2             !!!  Dipole parameterization  !!!
        dip2 = dip*dip

        xb = q2/(q2+w2-mp2)
        xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
        xpr = 1./xpr
c        t = log(log((q2+xval(50))/0.330**2)/log(xval(50)/0.330**2))

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
          mass(7) = xval(41)
          intwidth(7) = xval(42)
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

          pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)

          if(i.EQ.2) then
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif 

          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

          pgam(i) = intwidth(i)*pgam(i)

          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)

        enddo
 

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo

          if(sf.EQ.1) then

            if(i.eq.6) height(i) = rescoef(i,1)/
     &        (1.+ q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)


             height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.71)**rescoef(i,4)

          else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif

        enddo

CCC    End resonance Q^2 dependence calculations   CCC

     
        do i=1,3               !!!  Non-Res coefficients  !!!
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo


        if(sf.EQ.2) then      !!!  4th resonance region  !!!
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
        else
          height(7) = xval(49)*dip2 
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC

        sig_res = 0.0

        do i=1,7
          sigr(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(i) = sigr(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
          sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
          sig_res = sig_res + sigr(i)   
        enddo


CCC    Finish resonances / start non-res background calculation   CCC

 
        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          enddo

          sig_nr = sig_nr*xpr


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xpr)**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif


        sig = sig_res + sig_nr

          
        if(L.EQ.1) sigtemp = sig  

      enddo
       

CCC   Now extrapolate sig_L linearly to zero for Q^2 less than q2min   CCC

      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/dq2
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif


 1000  format(9f12.5)

      RETURN 
      END 

! Christy fit to proton
      SUBROUTINE CHRISTY806(W2,Q2,F1,R,sigt,sigl)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r,W1p,W2p,nu
      real*8 noverp,fnp_nmc
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval1 / 
     & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
     & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
     & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
     & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
     & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
     & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
     & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
     & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
     & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
     & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


      data xvalL/
     & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
     & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
     & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
     & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
     & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
     & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
     & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
     & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
     & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
     & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

      W1p=0.
      W2p=0.
      R=0.
      F1=0.
      sigl=0.
      sigt=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmod316(1,w2,q2,xval1,sigT)
      call resmod316(2,w2,q2,xvalL,sigL)


      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

      NU = q2 / 2. / mp / xb
      W1p = F1 / MP 
      W2p = W1p /(1.0 + NU*NU/Q2) * (1.0 + R)
      return
   
      end

CCC  Version 031606  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and fLparms.dat.  Units are ub/Sr/Gev.                          CCC

             
      SUBROUTINE RESMOD316(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low,dq2,t,xpr
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2
      common/tst1/sigr,sig_nr

      lowq2 = .false.
      lmax = 1
      q2temp = q2
      dq2 = 0.05

      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w

      br(1,1) = 1.0     !!! single pion branching ratios
      br(2,1) = 0.5
      br(3,1) = 0.65
      br(4,1) = 0.65
      br(5,1) = 0.4
      br(6,1) = 0.65
      br(7,1) = 0.6

      if(sf.EQ.2) then 
        br(6,1) = xval(48)
        br(2,1) = xval(49)
      endif 

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  ? 4th resonance region

      do i=1,7
        x0(i) = 0.165
      enddo
      x0(4) = 0.6

      do i=1,7
        br(i,2) = 1.-br(i,1)
      enddo
    

      if(sf.EQ.1) then
        q2low = 0.00
      else
        q2low = 0.1
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + dq2
        endif

        dip = 1./(1.+q2/0.71)**2             !!!  Dipole parameterization  !!!
        dip2 = dip*dip

        xb = q2/(q2+w2-mp2)
        xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
        xpr = 1./xpr
c        t = log(log((q2+xval(50))/0.330**2)/log(xval(50)/0.330**2))

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
          mass(7) = xval(41)
          intwidth(7) = xval(42)
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

          pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)

          if(i.EQ.2) then
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif 

          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

          pgam(i) = intwidth(i)*pgam(i)

          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)

        enddo
 

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo

          if(sf.EQ.1) then

            if(i.eq.6) height(i) = rescoef(i,1)/
     &        (1.+ q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)


             height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.71)**rescoef(i,4)

          else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif

        enddo

CCC    End resonance Q^2 dependence calculations   CCC

     
        do i=1,3               !!!  Non-Res coefficients  !!!
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo


        if(sf.EQ.2) then      !!!  4th resonance region  !!!
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
        else
          height(7) = xval(49)*dip2 
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC

        sig_res = 0.0

        do i=1,7
          sigr(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(i) = sigr(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
          sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
          sig_res = sig_res + sigr(i)   
        enddo


CCC    Finish resonances / start non-res background calculation   CCC

 
        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          enddo

          sig_nr = sig_nr*xpr


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xpr)**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif


        sig = sig_res + sig_nr

          
        if(L.EQ.1) sigtemp = sig  

      enddo
       

CCC   Now extrapolate sig_L linearly to zero for Q^2 less than q2min   CCC

      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/dq2
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif


 1000  format(9f12.5)

      RETURN 
      END 



      SUBROUTINE ARENHOVEL_INIT()
      IMPLICIT NONE
      REAL*8 QSQ_MAX,QSQ_GEV,NU,E_NP_FIND,Q2_CM_TEMP,QSQ_CM_MAX,WSQ,NU_T
      INTEGER NQ_,NE_,NMODEL
c      PARAMETER(NQ_=26,NE_=90,NMODEL=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL=2)
      INTEGER NE(0:NQ_,NMODEL)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL),Q2_CM(0:NE_,0:NQ_,NMODEL),
     >  FLT(2,0:NE_,0:NQ_,NMODEL),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE

c      CHARACTER IN_FILE(NMODEL)*100,QSQ_STR*4
      CHARACTER IN_FILE(NMODEL)*100,QSQ_STR*4,qsq_str2*2
      INTEGER IQ,IE,IM
      REAL*8 CHBAR,ALPHA,PI
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PI    = 3.1415927)   
      REAL*8 CONV/25.6819/ ! GeV^2 to fm_2
c      REAL*8 QSQ_V(0:NQ_)/0.,.125, .25, .5, .75, 1.0, 1.25, 1.5, 1.75,
c     >   2., 2.25, 2.5, 3.0, 3.5, 4.0,  4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
c     >   7.5, 8.0, 8.5, 9.0, 9.5, 10./
      REAL*8 QSQ_V(0:NQ_)/0.,.125, .25, .5, .75, 1.0, 1.25, 1.5, 1.75,
     >   2., 2.25, 2.5, 3.0, 3.5, 4.0,  4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
     >   7.5, 8.0, 8.5, 9.0, 9.5, 10., 11., 12., 13., 14., 15., 16.,
     >  17., 18., 19., 20., 22., 24., 26., 28., 30., 32., 34., 36.,
     >  38., 40. /
      REAL*8 E_NP_ADD/20./! extrapolate E_np this many MeV, to F=0
  

! cross section in nb/sr/GeV
      CMOTT  = 0.389E6 * ALPHA**2/4.  ! pb  
      DO IQ=1,NQ_
       QSQ_VAL(IQ) = QSQ_V(IQ) ! put into common block
       IF(QSQ_VAL(IQ).LT.0.2) THEN
        WRITE(QSQ_STR,'(F4.3)')QSQ_VAL(IQ)
       ELSEIF(QSQ_VAL(IQ).GE.10.0) THEN
        WRITE(QSQ_STR,'(F4.1)')QSQ_VAL(IQ)
       ELSE
        WRITE(QSQ_STR,'(F4.2)')QSQ_VAL(IQ)
       ENDIF
       IN_FILE(1)='input/cross-sec-data/arenhovel_v18_'
     >   //QSQ_STR//'_T.dat'
       IN_FILE(2) ='input/cross-sec-data/arenhovel_bonn_' 
     >   //QSQ_STR//'_T.dat'
       if(QSQ_VAL(IQ).GT.10.0) THEN
         WRITE(QSQ_STR2,'(i2.2)')int(QSQ_VAL(IQ))
         IN_FILE(1)='input/cross-sec-data/Bosted_v18_'
     >     //QSQ_STR2//'_T.dat'
         IN_FILE(2)='input/cross-sec-data/Bosted_bonn_'
     >     //QSQ_STR2//'_T.dat'
       endif
       DO IM=1,NMODEL
        OPEN(9,FILE=IN_FILE(IM))
        READ(9,'(/)') ! skip 2 header lines
        DO IE=1,NE_
         if(QSQ_VAL(IQ).GT.10.0) THEN
           E_NP(IE,IQ,IM) = E_NP(IE,25,IM)
           READ(9,'(27X,2E11.3)',END=99) FLT(1,IE,IQ,IM),FLT(2,IE,IQ,IM)
           NE(IQ,IM)=IE +1 ! extra 1 for extrapolation beyond maximum Enp
         else
          READ(9,'(2X,4E11.3)',END=99) E_NP(IE,IQ,IM), 
     >     Q2_CM(IE,IQ,IM),FLT(1,IE,IQ,IM),FLT(2,IE,IQ,IM)
          NE(IQ,IM)=IE +1 ! extra 1 for extrapolation beyond maximum Enp
         endif
        ENDDO ! ie
 99     close(unit=9)
       ENDDO ! im
      ENDDO ! iq
! Put in zero values for QSQ=0
      DO IM=1,NMODEL
       NE(0,IM) = NE(1,IM)     
       DO IE=0,NE(0,IM)
        FLT(1,IE,0,IM) = 0.
        FLT(2,IE,0,IM) = 0.
        E_NP(IE,0,IM) = E_NP(IE,1,IM)
       ENDDO ! ie

       DO IQ=1,NQ_
! Put in zero values for Enp=0
        FLT(1,0,IQ,IM) = 0.
        FLT(2,0,IQ,IM) = 0.
        E_NP(0,IQ,IM) =  0.
       ENDDO
       DO IQ=0,NQ_
! Put in zero values for Enp > Enp Max (beyond quasi-elastic peak)
        E_NP(NE(IQ,IM),IQ,IM) =  E_NP(NE(IQ,IM)-1,IQ,IM) + E_NP_ADD
        E_NP(NE(IQ,IM)+1,IQ,IM) = E_NP(NE(IQ,IM),IQ,IM) +  E_NP_ADD
        FLT(1,NE(IQ,IM),IQ,IM) = 0. 
        FLT(1,NE(IQ,IM)+1,IQ,IM) = 0.
        FLT(2,NE(IQ,IM),IQ,IM) = 0.
        FLT(2,NE(IQ,IM)+1,IQ,IM) = 0.

       ENDDO ! iq
      ENDDO ! im
      QSQ_CM_MAX = QSQ_VAL(NQ_)/CONV
      RETURN
      END

!================================================================
      real*8 function  ARENHOVEL_SIG(E,EP,TH,MODEL,IN_RANGE)
!-----------------------------------------------
! Get cross section from the structure functions
!-----------------------------------------------
      IMPLICIT NONE
      REAL*4 E,EP,TH
      INTEGER MODEL
      INTEGER IN_RANGE
      INTEGER NQ_,NE_,NMODEL_
c      PARAMETER(NQ_=26,NE_=90,NMODEL_=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL_=2)
      INTEGER NE(0:NQ_,NMODEL_)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL_),Q2_CM(0:NE_,0:NQ_,NMODEL_),
     >  FLT(2,0:NE_,0:NQ_,NMODEL_),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE
      REAL*8 QSQ,W1,W2,SIN2T2,CSMOTT,NU
     

      SIN2T2 =(SIN(.017453*TH/2))**2
      QSQ = 4.*E*EP*SIN2T2
      NU = E -EP

      CALL ARENHOVEL_W12(QSQ,NU,W1,W2,MODEL,IN_RANGE)

      CSMOTT =  CMOTT/(E*SIN2T2)**2 
      ARENHOVEL_SIG = CSMOTT *(W2*(1.-SIN2T2) +2.*W1 *SIN2T2)
      RETURN
      END

!=========================================================================
      SUBROUTINE ARENHOVEL_W12(QSQ,NU,W1,W2,MODEL,IN_RANGE)
!---------------------------------------------
! Get W1 and W2 from F_L and F_T
!--------------------------------------------
      IMPLICIT NONE     
      REAL*8 QSQ,NU,W1,W2 
      INTEGER MODEL
      INTEGER IN_RANGE
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/
      REAL*8 F_OUT(2),QSQ_L,Wnp2
      REAL*8 PI2SQ /19.7392/ ! 2*PI**2
      REAL*8 GEV_FM/5.07/ ! 1/(GeV-F 
      REAL*8 ALPHA
      PARAMETER (ALPHA = 7.29735E-03)    

      CALL ARENHOVEL_CALC(QSQ,NU,F_OUT,MODEL,IN_RANGE)
      QSQ_L = QSQ +NU**2  ! 3 momentum squared (lab)
      Wnp2 = MD*(MD+2.*NU) -QSQ
      W2 = (Wnp2/MDSQ*(QSQ/QSQ_L)**2 *F_OUT(1) +
     >      .5*QSQ/QSQ_L *F_OUT(2))/PI2SQ * GEV_FM/ALPHA
      W1 = .5*F_OUT(2)/PI2SQ * GEV_FM/ALPHA

      RETURN
      END

!==================================================================
!--------------------------------------------------------
! Get F_L [F_out(1)] and F_T [F_OUT(2)] by interpolating tables 
!  provided by Arenhovel.
! Lowest Q2 =.125 fm^-2 : Interpolate linearly to zero Q2 and zero F's
! Interpolate linearly from lowest Enp given to Enp=0 and zero F's
! ERROR codes in IN_RANGE
! = 0 OK
! = 1 E_NP < 0    sets output to 0 (correct)
! = 2 Enp > Max in Tables. extrapolate to zero at 20+Enp Max    (OK)
! = 3 QSQ > max in tables
! = 4 QSQ > max in tables AND E_NP > MAX in tables
! <3 is OK to use
! >=3 nead other fit for higher Q^2
!-----------------------------------------------------------------

      SUBROUTINE ARENHOVEL_CALC(QSQ,NU,F_OUT,MODEL,IN_RANGE)
! Get F_L [F_out(1)] and F_T [F_OUT(2)] by interpolating tables 
!-----------------------------------------------------------
      IMPLICIT NONE
      REAL*8 QSQ,NU,F_OUT(2) ! f_out(1) = F_L, f_out(2) = F_t
      INTEGER MODEL
      INTEGER IN_RANGE ! 0= OK, 1=Enp<0, 2= Enp high, 3=Q2 hi, 4=Q2>max && E_NP>MAX 
    
      INTEGER NQ_,NE_,NMODEL_,JQ1,JQ2
c      PARAMETER(NQ_=26,NE_=90,NMODEL_=2)
      PARAMETER(NQ_=46,NE_=90,NMODEL_=2)
      INTEGER NE(0:NQ_,NMODEL_)
      REAL*8 E_NP(0:NE_,0:NQ_,NMODEL_),Q2_CM(0:NE_,0:NQ_,NMODEL_),
     >  FLT(2,0:NE_,0:NQ_,NMODEL_),QSQ_VAL(0:NQ_),CMOTT
      COMMON/ARENHOVEL/E_NP,Q2_CM,FLT,QSQ_VAL,CMOTT,NE

      INTEGER IQ,JQ,je(2),KE,LT
      REAL*8 Q2_CM_GEV,E_NP_FIND,Q2_CM_FIND
      REAL*8  denomq,denome(2),w1q,w2q,w1e(2),w2e(2),fac1,fac2,w1ee,w2ee
      REAL*8 CONV/25.6819/ ! GeV^2 to fm_2

      IN_RANGE =0 

      CALL  KIN_TRANSFORM_NORM_ARENHOVEL(QSQ,NU,E_NP_FIND,Q2_CM_GEV)
      IF(E_NP_FIND.LT.0.) THEN
       IN_RANGE =1
       GOTO 999
      ENDIF
      Q2_CM_FIND = Q2_CM_GEV *CONV  ! convert to fm^{-2}
      E_NP_FIND = E_NP_FIND *1000.   ! convert to MeV from GeV

      IF(Q2_CM_FIND.GT.QSQ_VAL(NQ_))THEN ! QSQ > max in tables
       IN_RANGE=3
       IF(E_NP_FIND.GT.E_NP(NE(NQ_,MODEL),NQ_,MODEL)) THEN
        IN_RANGE=4
        GO TO 999
       ENDIF
       IQ = NQ_-1  ! QSQ > max in tables
       W1Q =1.
       W2Q =0.
       DENOMQ=1.
       JQ1=IQ
       JQ2=IQ
       GOTO 1100
      ENDIF
      DO iq=0, NQ_-1
       if((Q2_CM_FIND.GE.QSQ_VAL(iq)).AND.
     >    (Q2_CM_FIND.LE.QSQ_VAL(iq+1)))      THEN
        denomq =ABS(QSQ_VAL(iq+1) -QSQ_VAL(iq))
        w1q = abs(Q2_CM_FIND-QSQ_VAL(iq))
        w2q = abs(QSQ_VAL(iq+1)-Q2_CM_FIND)
        GO TO 1000
       ENDIF
      ENDDO  ! iq
1000  CONTINUE
      JQ1=IQ
      JQ2=IQ+1
      IF((E_NP_FIND.GT.E_NP(NE(IQ,MODEL),IQ,MODEL)).AND.
     >   (E_NP_FIND.GT.E_NP(NE(IQ+1,MODEL),IQ+1,MODEL))) THEN
       IN_RANGE=2 
       GOTO 999
      ENDIF
 1100 CONTINUE
! find E_NP factors for each of the flanking QSQ's
      DO JQ=JQ1,JQ2
       if(jq.gt.nq_ .or.jq.lt.0) then
         write(6,'(''ERROR!'',4i4,4f10.3)') iq,jq,jq1,jq2,
     >     Q2_CM_FIND.GE.QSQ_VAL(nq_),QSQ_VAL(0),E_NP_FIND
       endif
       if(E_NP_FIND.LE.E_NP(0,JQ,MODEL)) THEN !e_np(0,jq,model) =0 always
        w1e(JQ-IQ+1)=0.
        w2e(JQ-IQ+1)=1.
        denome(JQ-IQ+1)=1.
        je(jq-iq+1)=0
       elseif(E_NP_FIND.GE.E_NP(NE(JQ,MODEL),JQ,MODEL)) THEN !Enp > max
        w1e(JQ-IQ+1)=1.
        w2e(JQ-IQ+1)=0.
        denome(JQ-IQ+1)=1.
        je(jq-iq+1)=NE(JQ,MODEL)
       else 
        DO ke=0, NE(JQ,MODEL) -1 
         if((E_NP_FIND.GE.E_NP(ke,jq,MODEL)).AND.
     >     (E_NP_FIND.LE.E_NP(ke+1,jq,MODEL)))      THEN
          denome(JQ-IQ+1) =
     >     ABS(E_NP(ke+1,jq,MODEL) -E_NP(ke,jq,MODEL))
          w1e(JQ-IQ+1) = abs(E_NP_FIND -E_NP(ke,jq,MODEL))
          w2e(JQ-IQ+1) = abs(E_NP(ke+1,jq,MODEL)-E_NP_FIND)
          je(JQ-IQ+1) = ke
          go to 2000
         endif
        ENDDO ! ke
       ENDIF !e_np_find
 2000  continue
      ENDDO  ! jq
      IF(JQ1.EQ.JQ2) THEN  ! high edge of Enp_find range
        JE(2) = JE(1)
        W1E(2) = W1E(1)
        W2E(2) = W2E(1)
        DENOME(2) = DENOME(1)
      ENDIF
      w1ee = (w2q*w1e(1)/denome(1) + w1q*w1e(2)/denome(2))/denomq
      w2ee = (w2q*w2e(1)/denome(1) + w1q*w2e(2)/denome(2))/denomq
          
      DO LT =1,2
       fac1 = 
     >  (w1q*FLT(lt,je(2),iq+1,MODEL)+w2q*FLT(lt,je(1),iq,MODEL))/denomq
       fac2 = 
     >  (w1q*FLT(lt,je(2)+1,iq+1,MODEL)+w2q*FLT(lt,je(1)+1,iq,MODEL))/
     >   denomq
       F_OUT(LT) =w1ee*fac2 + w2ee*fac1
      ENDDO
      RETURN

 999  CONTINUE
      F_OUT(1)=0.
      F_OUT(2)=0.
      RETURN
      END

      SUBROUTINE KIN_TRANSFORM_ARENHOVEL_NORM (E_NP,Q2_CM,WSQ,QSQ,NU)
!------------------------------------------------------
! transform form Arenhovel's coordinates of E_NP, q^2_cm to 
! what I normally use.
!-------------------------------------------------------
      IMPLICIT NONE
      REAL*8 E_NP, Q2_CM, WSQ, QSQ,NU
      REAL*8 C,ARG
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/

      
      C = MDSQ  - (E_NP +MN+MP)**2  -Q2_CM*(E_NP +MN+MP)**2/MDSQ
      ARG = 4.*MDSQ -4.*C
      NU = (-MD2 + SQRT(ARG))/2.
      QSQ = MDSQ+MD2*NU -(E_NP+MN+MP)**2
      WSQ = MPSQ +2.*MP*NU -QSQ
      RETURN
      END
    
      SUBROUTINE KIN_TRANSFORM_NORM_ARENHOVEL(QSQ,NU,E_NP,Q2_CM)   
!------------------------------------------------------
! Everything in GeV
! transform form what I normally use.QSQ,NU
! to Arenhovel's coordinates of E_NP, q^2_cm
!-------------------------------------------------------
      IMPLICIT NONE
      REAL*8 E_NP, Q2_CM, QSQ,NU
      REAL*8 MD/1.87561/
      REAL*8 MDSQ/3.51788/
      REAL*8 MD2/3.75122/
      REAL*8 MP/.93827/,MN/.93957/
      REAL*8 MPSQ/.88035/
      REAL*8 Wnp2

      Wnp2  = MDSQ +MD2*NU -QSQ                 ! eq 5 Arenhovel notes
      E_NP = SQRT(Wnp2) -MP -MN   ! eq. 1. Arenhovel notes
      Q2_CM = MDSQ*(NU**2+QSQ)/Wnp2
      RETURN
      END
